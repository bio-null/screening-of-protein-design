import os
import argparse
import mdtraj as md
import logging
import shutil

def setup_logging(root_path, log_level=logging.INFO):
    os.makedirs(root_path, exist_ok=True)
    log_file = os.path.join(root_path, "global_rmsd_filter.log")
    
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger("global_rmsd_filter")

def global_rmsd(pdb1, pdb2):
    try:
        trajectory1 = md.load(pdb1)
        trajectory2 = md.load(pdb2)
        
        atoms1 = trajectory1.topology.select("name CA")
        atoms2 = trajectory2.topology.select("name CA")
        
        if len(atoms1) != len(atoms2):
            raise ValueError(f"Cα原子数量不同: {len(atoms1)} vs {len(atoms2)}")
        
        trajectory1.superpose(trajectory2, atom_indices=atoms1, ref_atom_indices=atoms2)
        
        rmsd_value = md.rmsd(trajectory1, trajectory2, atom_indices=atoms1, ref_atom_indices=atoms2)
        
        return rmsd_value[0]
    
    except Exception as e:
        logging.error(f"计算 {pdb1} 和 {pdb2} 的RMSD时出错: {str(e)}")
        return None

def filter_pdbs_by_rmsd(design_path, reference_pdb, output_dir, rmsd_threshold, logger):
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_files = [f for f in os.listdir(design_path) if f.endswith(".pdb")]
    total_files = len(pdb_files)
    passed_files = 0
    
    logger.info(f"开始处理目录: {design_path}")
    logger.info(f"找到 {total_files} 个PDB文件")
    logger.info(f"参考结构: {reference_pdb}")
    logger.info(f"RMSD阈值: {rmsd_threshold} nm")
    
    for pdb_file in pdb_files:
        file_path = os.path.join(design_path, pdb_file)
        
        rmsd_value = global_rmsd(file_path, reference_pdb)
        
        if rmsd_value is None:
            logger.warning(f"无法计算 {pdb_file} 的RMSD，跳过此文件")
            continue
            
        rmsd_angstrom = rmsd_value * 10.0
        
        logger.info(f"文件: {pdb_file} - RMSD: {rmsd_value:.3f} nm ({rmsd_angstrom:.2f} Å)")
        
        if rmsd_value <= rmsd_threshold:
            dest_path = os.path.join(output_dir, pdb_file)
            shutil.copy(file_path, dest_path)
            logger.info(f"  符合条件! 已复制到: {dest_path}")
            passed_files += 1
        else:
            logger.info(f"  不符合条件 (要求 <= {rmsd_threshold:.3f} nm)")
    
    logger.info("\n===== 处理完成 =====")
    logger.info(f"总文件数: {total_files}")
    logger.info(f"符合条件文件数: {passed_files}")
    
    if total_files > 0:
        logger.info(f"通过率: {passed_files/total_files*100:.2f}%")
    else:
        logger.warning("未找到任何PDB文件")
    
    return passed_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据与参考结构的RMSD差异筛选蛋白质设计结果")
    parser.add_argument(
        "-d", "--design-path",
        type=str,
        required=True,
        help="包含设计PDB文件的文件夹路径"
    )
    parser.add_argument(
        "-r", "--reference-pdb",
        type=str,
        required=True,
        help="参考PDB文件路径（原始序列结构）"
    )
    parser.add_argument(
        "-t", "--rmsd-threshold",
        type=float,
        default=0.2,
        help="RMSD过滤阈值（单位：纳米，默认0.2nm=2Å）"
    )
    parser.add_argument(
        "--root-path",
        default="/home/ug2023/ug523111910012/project1/filter/rmsd",
        help="根目录路径"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        help="输出目录路径（可选，默认在根路径下创建）"
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="日志详细程度"
    )
    
    args = parser.parse_args()
    
    design_path = os.path.normpath(args.design_path)
    reference_pdb = os.path.normpath(args.reference_pdb)
    root_path = os.path.normpath(args.root_path)
    rmsd_threshold = args.rmsd_threshold
    
    log_level = getattr(logging, args.log_level)
    logger = setup_logging(root_path, log_level)
    
    if args.output_dir:
        output_dir = os.path.normpath(args.output_dir)
    else:
        output_dir = os.path.join(root_path, f"filtered_rmsd_{rmsd_threshold:.2f}".replace('.', '_'))
    
    logger.info(f"输出目录设置为: {output_dir}")
    
    if not os.path.exists(design_path):
        logger.error(f"设计路径不存在: {design_path}")
        exit(1)
    
    if not os.path.isdir(design_path):
        logger.error(f"设计路径不是目录: {design_path}")
        exit(1)
    
    if not os.path.exists(reference_pdb):
        logger.error(f"参考PDB文件不存在: {reference_pdb}")
        exit(1)
    
    logger.info(f"开始筛选过程...")
    logger.info(f"设计路径: {design_path}")
    logger.info(f"参考结构: {reference_pdb}")
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"RMSD阈值: {rmsd_threshold} nm ({rmsd_threshold*10:.1f} Å)")
    
    passed_count = filter_pdbs_by_rmsd(design_path, reference_pdb, output_dir, rmsd_threshold, logger)
    
    logger.info(f"筛选完成! 符合条件的文件数: {passed_count}")
    logger.info(f"输出目录: {output_dir}")