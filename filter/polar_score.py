import os
import argparse
import logging
import shutil
import math
import torch
import mdtraj as md

def setup_logging(root_path, log_level=logging.INFO):
    os.makedirs(root_path, exist_ok=True)
    log_file = os.path.join(root_path, "surface_polar_filter.log")
    
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger("surface_polar_filter")

def load_pdb(pdb_file):
    traj = md.load(pdb_file)
    return traj.topology, traj

def calculate_backbone_neighbors(pdb_file):
    top, traj = load_pdb(pdb_file)
    CBs_index = top.select("backbone and name C")
    CBs_coords = traj.xyz[0, CBs_index]
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    CBs_coords = torch.from_numpy(CBs_coords).to(device)
    CBs_dis = torch.norm(CBs_coords[:, None] - CBs_coords, dim=2, p=2)
    return CBs_dis

def classify_polarity(traj):
    residues = [str(residue) for residue in traj.topology.residues]
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    non_polar = torch.zeros(len(residues), device=device)
    polar = torch.zeros(len(residues), device=device)
    
    non_polar_residues = ["ILE", "LEU", "MET", "TRP", "PHE", "VAL"]
    polar_residues = ["SER", "THR", "TYR", "ASN", "GLN"]
    
    for i, residue in enumerate(residues):
        res_name = residue[:3]
        if res_name in non_polar_residues:
            non_polar[i] = 1
        elif res_name in polar_residues:
            polar[i] = 1
            
    return non_polar, polar

def calculate_angles(traj):
    top = traj.topology
    CBs_index = top.select("backbone and name C")
    CBs_coords = traj.xyz[0, CBs_index]
    CAs_index = top.select("backbone and name CA")
    CAs_coords = traj.xyz[0, CAs_index]
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    CAs_CBs_coords = torch.from_numpy(CBs_coords - CAs_coords).to(device)
    
    norm = torch.norm(CAs_CBs_coords, dim=1, keepdim=True)
    normalized = CAs_CBs_coords / norm
    
    cos_angles = torch.mm(normalized, normalized.t())
    cos_angles = torch.clamp(cos_angles, -1.0, 1.0)
    phi_ij = torch.acos(cos_angles)
    
    return phi_ij

def calculate_surface_polar_score(pdb_file):
    try:
        top, traj = load_pdb(pdb_file)
        
        CBs_dis = calculate_backbone_neighbors(pdb_file)
        
        non_polar, polar = classify_polarity(traj)
        
        phi_ij = calculate_angles(traj)
        
        m = 1.0
        a = 0.5
        b = 2.0
        
        distance_weights = 1 / (1 + torch.exp(CBs_dis - m))
        
        angle_weights = ((torch.cos(math.pi - phi_ij) + a) / (1 + a)) ** b
        
        combined_weights = distance_weights * angle_weights
        diag_mask = torch.eye(combined_weights.size(0), device=combined_weights.device).bool()
        combined_weights = combined_weights.masked_fill(diag_mask, 0)
        
        n_i = torch.sum(combined_weights, dim=1)
        
        median_n_i = torch.median(n_i)
        sigmoid_term = 1 - torch.sigmoid(n_i - median_n_i)
        
        numerator = torch.sum(non_polar * sigmoid_term)
        denominator = torch.sum(sigmoid_term)
        
        polar_score = numerator / denominator
        
        return polar_score.item()
    
    except Exception as e:
        logging.error(f"计算 {pdb_file} 表面极性分数时出错: {str(e)}")
        return None

def filter_pdbs_by_polar_score(
    design_path, 
    output_dir, 
    polar_threshold, 
    logger
):
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_files = [f for f in os.listdir(design_path) if f.endswith(".pdb")]
    total_files = len(pdb_files)
    passed_files = 0
    
    logger.info(f"开始处理目录: {design_path}")
    logger.info(f"找到 {total_files} 个PDB文件")
    logger.info(f"表面极性分数阈值: {polar_threshold:.4f}")
    
    for pdb_file in pdb_files:
        file_path = os.path.join(design_path, pdb_file)
        
        polar_score = calculate_surface_polar_score(file_path)
        
        if polar_score is None:
            logger.warning(f"无法计算 {pdb_file} 的表面极性分数，跳过此文件")
            continue
            
        logger.info(f"文件: {pdb_file} - 表面极性分数: {polar_score:.4f}")
        
        if polar_score <= polar_threshold:
            dest_path = os.path.join(output_dir, pdb_file)
            shutil.copy(file_path, dest_path)
            logger.info(f"  符合条件! 已复制到: {dest_path}")
            passed_files += 1
        else:
            logger.info(f"  不符合条件 (要求 <= {polar_threshold:.4f})")
    
    logger.info("\n===== 处理完成 =====")
    logger.info(f"总文件数: {total_files}")
    logger.info(f"符合条件文件数: {passed_files}")
    
    if total_files > 0:
        logger.info(f"通过率: {passed_files/total_files*100:.2f}%")
    else:
        logger.warning("未找到任何PDB文件")
    
    return passed_files

def calculate_single_pdb(pdb_file, logger):
    if not os.path.exists(pdb_file):
        logger.error(f"PDB文件不存在: {pdb_file}")
        return
    
    polar_score = calculate_surface_polar_score(pdb_file)
    
    if polar_score is not None:
        logger.info(f"文件: {os.path.basename(pdb_file)} - 表面极性分数: {polar_score:.4f}")
        print(f"表面极性分数: {polar_score:.4f}")
    else:
        logger.error(f"无法计算 {pdb_file} 的表面极性分数")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--design-path",
        type=str,
        help="包含设计PDB文件的文件夹路径"
    )
    parser.add_argument(
        "-r", "--reference-pdb",
        type=str,
        help="参考PDB文件路径"
    )
    parser.add_argument(
        "-p", "--polar-score",
        type=float,
        help="表面极性分数过滤阈值"
    )
    parser.add_argument(
        "-c", "--calculate",
        type=str,
        help="计算单个PDB文件的表面极性分数"
    )
    parser.add_argument(
        "--root-path",
        type=str,
        default="/home/ug2023/ug523111910012/project1/filter/polar_score",
        help="根目录路径"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        help="输出目录路径（可选）"
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="日志详细程度"
    )
    
    args = parser.parse_args()
    
    if args.calculate:
        mode = "single"
    elif args.design_path:
        if args.reference_pdb and args.polar_score:
            parser.error("不能同时使用 -r 和 -p 参数")
        elif not args.reference_pdb and not args.polar_score:
            parser.error("必须提供阈值")
        else:
            mode = "filter"
    else:
        parser.error("必须提供 -d 或 -c 参数")
    
    root_path = args.root_path
    log_level = getattr(logging, args.log_level)
    logger = setup_logging(root_path, log_level)
    
    if mode == "single":
        calculate_single_pdb(args.calculate, logger)
    else:
        design_path = args.design_path
        
        if args.reference_pdb:
            reference_pdb = args.reference_pdb
            polar_threshold = calculate_surface_polar_score(reference_pdb)
            if polar_threshold is None:
                logger.error(f"无法计算参考PDB {reference_pdb} 的表面极性分数")
                exit(1)
            logger.info(f"参考PDB {reference_pdb} 的表面极性分数: {polar_threshold:.4f}")
        else:
            polar_threshold = args.polar_score
        
        if args.output_dir:
            output_dir = args.output_dir
        else:
            threshold_str = f"{polar_threshold:.4f}".replace('.', '_')
            output_dir = os.path.join(root_path, f"filtered_polar_{threshold_str}")
        
        logger.info(f"输出目录: {output_dir}")
        
        if not os.path.exists(design_path):
            logger.error(f"设计路径不存在: {design_path}")
            exit(1)
        
        if not os.path.isdir(design_path):
            logger.error(f"设计路径不是目录: {design_path}")
            exit(1)
        
        logger.info(f"开始筛选...")
        logger.info(f"设计路径: {design_path}")
        logger.info(f"输出目录: {output_dir}")
        logger.info(f"表面极性分数阈值: {polar_threshold:.4f}")
        
        passed_count = filter_pdbs_by_polar_score(
            design_path, output_dir, polar_threshold, logger
        )
        
        logger.info(f"筛选完成! 符合条件的文件数: {passed_count}")
        logger.info(f"输出目录: {output_dir}")