import os
import argparse
import mdtraj as md
import logging
import shutil

def setup_logging(root_path, log_level=logging.INFO):
    os.makedirs(root_path, exist_ok=True)
    log_file = os.path.join(root_path, "sasa_filter.log")
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger("sasa_filter")

def compute_sasa(pdb_file):
    try:
        traj = md.load(pdb_file)
        sasa = md.shrake_rupley(traj)[0]
        return sasa.sum()
    except Exception as e:
        logging.error(f"计算SASA时出错: {str(e)}")
        return None

def filter_pdbs(design_path, output_dir, min_sasa, max_sasa, logger):
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = [f for f in os.listdir(design_path) if f.endswith(".pdb")]
    total_files = len(pdb_files)
    passed_files = 0
    
    logger.info(f"开始处理目录: {design_path}")
    logger.info(f"找到 {total_files} 个PDB文件")
    logger.info(f"SASA范围: {min_sasa:.1f} - {max_sasa:.1f} nm²")
    
    for pdb_file in pdb_files:
        file_path = os.path.join(design_path, pdb_file)
        sasa_value = compute_sasa(file_path)
        
        if sasa_value is None:
            logger.warning(f"无法计算 {pdb_file} 的SASA，跳过此文件")
            continue
            
        logger.info(f"文件: {pdb_file} - SASA: {sasa_value:.1f} nm²")
        
        if min_sasa <= sasa_value <= max_sasa:
            dest_path = os.path.join(output_dir, pdb_file)
            shutil.copy(file_path, dest_path)
            logger.info(f"  符合条件! 已复制到: {dest_path}")
            passed_files += 1
        else:
            logger.info(f"  不符合条件 (要求范围: {min_sasa:.1f} - {max_sasa:.1f} nm²)")
    
    logger.info("\n===== 处理完成 =====")
    logger.info(f"总文件数: {total_files}")
    logger.info(f"符合条件文件数: {passed_files}")
    if total_files > 0: logger.info(f"通过率: {passed_files/total_files*100:.2f}%")
    else: logger.warning("未找到任何PDB文件")
    return passed_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据溶剂可及表面积筛选蛋白质设计结果")
    parser.add_argument("-d", "--design-path", type=str, required=True, help="设计PDB文件夹路径")
    parser.add_argument("--min-sasa", type=float, default=100.0, help="最小SASA(nm²)")
    parser.add_argument("--max-sasa", type=float, default=200.0, help="最大SASA(nm²)")
    parser.add_argument("--root-path", default="./filter_results/sasa", help="根目录路径")
    parser.add_argument("--output-dir", type=str, help="自定义输出目录")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO", help="日志级别")
    
    args = parser.parse_args()
    design_path = os.path.normpath(args.design_path)
    root_path = os.path.normpath(args.root_path)
    min_sasa = args.min_sasa
    max_sasa = args.max_sasa
    log_level = getattr(logging, args.log_level)
    logger = setup_logging(root_path, log_level)
    
    output_dir = args.output_dir if args.output_dir else os.path.join(root_path, f"filtered_sasa_{min_sasa:.0f}to{max_sasa:.0f}")
    logger.info(f"输出目录: {output_dir}")
    
    if not os.path.exists(design_path):
        logger.error(f"设计路径不存在: {design_path}")
        exit(1)
    if not os.path.isdir(design_path):
        logger.error(f"设计路径不是目录: {design_path}")
        exit(1)
    
    logger.info(f"开始筛选...")
    passed_count = filter_pdbs(design_path, output_dir, min_sasa, max_sasa, logger)
    logger.info(f"筛选完成! 符合条件的文件数: {passed_count}")