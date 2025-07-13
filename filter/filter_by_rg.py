import os
import argparse
import mdtraj as md
import logging
import shutil

def setup_logging(root_path, log_level=logging.INFO):
    os.makedirs(root_path, exist_ok=True)
    log_file = os.path.join(root_path, "rg_filter.log")
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger("rg_filter")

def radius_of_gyration(pdb_file):
    try:
        traj = md.load(pdb_file)
        rg = md.compute_rg(traj)[0]
        return rg
    except Exception as e:
        logging.error(f"计算回旋半径时出错: {str(e)}")
        return None

def filter_pdbs(design_path, output_dir, min_rg, max_rg, logger):
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = [f for f in os.listdir(design_path) if f.endswith(".pdb")]
    total_files = len(pdb_files)
    passed_files = 0
    
    logger.info(f"开始处理目录: {design_path}")
    logger.info(f"找到 {total_files} 个PDB文件")
    logger.info(f"回旋半径范围: {min_rg:.2f} - {max_rg:.2f} nm")
    
    for pdb_file in pdb_files:
        file_path = os.path.join(design_path, pdb_file)
        rg_value = radius_of_gyration(file_path)
        
        if rg_value is None:
            logger.warning(f"无法计算 {pdb_file} 的回旋半径，跳过此文件")
            continue
            
        logger.info(f"文件: {pdb_file} - 回旋半径: {rg_value:.3f} nm")
        
        if min_rg <= rg_value <= max_rg:
            dest_path = os.path.join(output_dir, pdb_file)
            shutil.copy(file_path, dest_path)
            logger.info(f"  符合条件! 已复制到: {dest_path}")
            passed_files += 1
        else:
            logger.info(f"  不符合条件 (要求范围: {min_rg:.3f} - {max_rg:.3f} nm)")
    
    logger.info("\n===== 处理完成 =====")
    logger.info(f"总文件数: {total_files}")
    logger.info(f"符合条件文件数: {passed_files}")
    if total_files > 0: logger.info(f"通过率: {passed_files/total_files*100:.2f}%")
    else: logger.warning("未找到任何PDB文件")
    return passed_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据回旋半径筛选蛋白质设计结果")
    parser.add_argument("-d", "--design-path", type=str, required=True, help="设计PDB文件夹路径")
    parser.add_argument("--min-rg", type=float, default=1.0, help="最小回旋半径(nm)")
    parser.add_argument("--max-rg", type=float, default=2.0, help="最大回旋半径(nm)")
    parser.add_argument("--root-path", default="./filter_results/radius_of_gyration", help="根目录路径")
    parser.add_argument("--output-dir", type=str, help="自定义输出目录")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO", help="日志级别")
    
    args = parser.parse_args()
    design_path = os.path.normpath(args.design_path)
    root_path = os.path.normpath(args.root_path)
    min_rg = args.min_rg
    max_rg = args.max_rg
    log_level = getattr(logging, args.log_level)
    logger = setup_logging(root_path, log_level)
    
    output_dir = args.output_dir if args.output_dir else os.path.join(root_path, f"filtered_rg_{min_rg:.1f}to{max_rg:.1f}")
    logger.info(f"输出目录: {output_dir}")
    
    if not os.path.exists(design_path):
        logger.error(f"设计路径不存在: {design_path}")
        exit(1)
    if not os.path.isdir(design_path):
        logger.error(f"设计路径不是目录: {design_path}")
        exit(1)
    
    logger.info(f"开始筛选...")
    passed_count = filter_pdbs(design_path, output_dir, min_rg, max_rg, logger)
    logger.info(f"筛选完成! 符合条件的文件数: {passed_count}")