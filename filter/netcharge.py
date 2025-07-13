import os
import argparse
import mdtraj as md
import logging
import shutil

def setup_logging(root_path, log_level=logging.INFO):
    os.makedirs(root_path, exist_ok=True)
    log_file = os.path.join(root_path, "net_charge_filter.log")
    
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger("net_charge_filter")

def netcharge(pdb_file):
    try:
        trajectory = md.load(pdb_file)
        top = trajectory.topology
        
        arg, lys, asp, glu = 0, 0, 0, 0
        name_list = ["ASP", "LYS", "GLU", "ARG"]
        
        for residue in top.residues:
            name = residue.name[:3]
            
            if name == name_list[0]:
                asp += 1
            elif name == name_list[1]:
                lys += 1
            elif name == name_list[2]:
                glu += 1
            elif name == name_list[3]:
                arg += 1
        
        net_charge = (arg + lys) - (asp + glu)
        return net_charge
    
    except Exception as e:
        logging.error(f"计算 {pdb_file} 净电荷时出错: {str(e)}")
        return None

def filter_pdbs_by_charge(design_path, output_dir, nc_threshold):
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_files = [f for f in os.listdir(design_path) if f.endswith(".pdb")]
    total_files = len(pdb_files)
    passed_files = 0
    
    logging.info(f"开始处理目录: {design_path}")
    logging.info(f"找到 {total_files} 个PDB文件")
    logging.info(f"净电荷阈值: {nc_threshold}")
    
    for pdb_file in pdb_files:
        file_path = os.path.join(design_path, pdb_file)
        charge = netcharge(file_path)
        
        if charge is None:
            continue
            
        logging.info(f"文件: {pdb_file} - 净电荷: {charge}")
        
        if charge <= nc_threshold:
            dest_path = os.path.join(output_dir, pdb_file)
            shutil.copy(file_path, dest_path)
            logging.info(f"  符合条件! 已复制到: {dest_path}")
            passed_files += 1
        else:
            logging.info(f"  不符合条件 (要求 <= {nc_threshold})")
    
    logging.info("\n===== 处理完成 =====")
    logging.info(f"总文件数: {total_files}")
    logging.info(f"符合条件文件数: {passed_files}")
    
    if total_files > 0:
        logging.info(f"通过率: {passed_files/total_files*100:.2f}%")
    else:
        logging.warning("未找到任何PDB文件")
    
    return passed_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据净电荷过滤蛋白质设计结果")
    parser.add_argument(
        "-d", "--design_path",
        type=str,
        required=True,
        help="包含设计PDB文件的文件夹路径"
    )
    parser.add_argument(
        "-n", "--net-charge",
        type=float,
        default=-1,
        help="净电荷过滤阈值"
    )
    parser.add_argument(
        "--root-path",
        default="/home/ug2023/ug523111910012/project1/filter/net_charge",
        help="根目录路径"
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="日志详细程度"
    )
    
    args = parser.parse_args()
    
    design_path = os.path.normpath(args.design_path)
    root_path = os.path.normpath(args.root_path)
    nc_threshold = args.net_charge
    
    log_level = getattr(logging, args.log_level)
    logger = setup_logging(root_path, log_level)
    
    output_dir = os.path.join(root_path, f"filtered_charge_{nc_threshold}")
    
    if not os.path.exists(design_path):
        logger.error(f"设计路径不存在: {design_path}")
        exit(1)
    
    if not os.path.isdir(design_path):
        logger.error(f"设计路径不是目录: {design_path}")
        exit(1)
    
    logger.info(f"开始过滤过程...")
    logger.info(f"设计路径: {design_path}")
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"净电荷阈值: {nc_threshold}")
    
    passed_count = filter_pdbs_by_charge(design_path, output_dir, nc_threshold)
    
    logger.info(f"过滤完成! 符合条件的文件数: {passed_count}")
    logger.info(f"输出目录: {output_dir}")