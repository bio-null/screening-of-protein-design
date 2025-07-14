import os
import argparse
import mdtraj as md
import logging
import shutil

def setup_logging(root_path, log_level=logging.INFO):
    os.makedirs(root_path, exist_ok=True)
    log_file = os.path.join(root_path, "local_rmsd_filter.log")
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger("local_rmsd_filter")

def local_rmsd(pdb1, pdb2, selection):
    try:
        traj1 = md.load(pdb1)
        traj2 = md.load(pdb2)
        atoms1 = traj1.topology.select(selection)
        atoms2 = traj2.topology.select(selection)
        
        if len(atoms1) != len(atoms2):
            raise ValueError(f"选择区域原子数量不同: {len(atoms1)} vs {len(atoms2)}")
        
        rmsd_value = md.rmsd(traj1, traj2, atom_indices=atoms1, ref_atom_indices=atoms2)
        return rmsd_value[0]
    except Exception as e:
        logging.error(f"计算局部RMSD时出错: {str(e)}")
        return None

def filter_pdbs(design_path, reference_pdb, output_dir, rmsd_threshold, selection, logger):
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = [f for f in os.listdir(design_path) if f.endswith(".pdb")]
    total_files = len(pdb_files)
    passed_files = 0
    
    logger.info(f"开始处理目录: {design_path}")
    logger.info(f"找到 {total_files} 个PDB文件")
    logger.info(f"参考结构: {reference_pdb}")
    logger.info(f"局部区域选择: {selection}")
    logger.info(f"RMSD阈值: {rmsd_threshold} nm")
    
    for pdb_file in pdb_files:
        file_path = os.path.join(design_path, pdb_file)
        rmsd_value = local_rmsd(file_path, reference_pdb, selection)
        
        if rmsd_value is None:
            logger.warning(f"无法计算 {pdb_file} 的局部RMSD，跳过此文件")
            continue
            
        rmsd_angstrom = rmsd_value * 10.0
        logger.info(f"文件: {pdb_file} - 局部RMSD: {rmsd_value:.3f} nm ({rmsd_angstrom:.2f} Å)")
        
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
    if total_files > 0: logger.info(f"通过率: {passed_files/total_files*100:.2f}%")
    else: logger.warning("未找到任何PDB文件")
    return passed_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据局部RMSD筛选蛋白质设计结果")
    parser.add_argument("-d", "--design-path", type=str, required=True, help="设计PDB文件夹路径")
    parser.add_argument("-r", "--reference-pdb", type=str, required=True, help="参考PDB文件路径")
    parser.add_argument("-t", "--rmsd-threshold", type=float, default=0.1, help="RMSD阈值(单位: nm)")
    parser.add_argument("-s", "--selection", type=str, default="name CA", help="原子选择表达式")
    parser.add_argument("--root-path", default="./filter_results/local_rmsd", help="根目录路径")
    parser.add_argument("--output-dir", type=str, help="自定义输出目录")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO", help="日志级别")
    
    args = parser.parse_args()
    design_path = os.path.normpath(args.design_path)
    reference_pdb = os.path.normpath(args.reference_pdb)
    root_path = os.path.normpath(args.root_path)
    rmsd_threshold = args.rmsd_threshold
    selection = args.selection
    log_level = getattr(logging, args.log_level)
    logger = setup_logging(root_path, log_level)
    
    output_dir = args.output_dir if args.output_dir else os.path.join(root_path, f"filtered_localrmsd_{rmsd_threshold:.2f}".replace('.', '_'))
    logger.info(f"输出目录: {output_dir}")
    
    if not os.path.exists(design_path):
        logger.error(f"设计路径不存在: {design_path}")
        exit(1)
    if not os.path.isdir(design_path):
        logger.error(f"设计路径不是目录: {design_path}")
        exit(1)
    if not os.path.exists(reference_pdb):
        logger.error(f"参考PDB文件不存在: {reference_pdb}")
        exit(1)
    
    logger.info(f"开始筛选...")
    passed_count = filter_pdbs(design_path, reference_pdb, output_dir, rmsd_threshold, selection, logger)
    logger.info(f"筛选完成! 符合条件的文件数: {passed_count}")
