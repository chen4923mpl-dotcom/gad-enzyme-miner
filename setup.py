#!/usr/bin/env python3
"""
GAD酶挖掘工具的安装脚本
"""

import os
import sys

def install():
    """安装工具"""
    print("GAD酶基因组挖掘工具安装")
    print("==========================")
    
    # 检查Python版本
    if sys.version_info < (3, 7):
        print("需要Python 3.7或更高版本")
        return False
    
    # 检查是否已有文件
    files = [
        "optimized_gad_mine.py",
        "gad_mine_simplified.py",
        "gad_config.json",
        "README.md",
        "requirements.txt"
    ]
    
    missing_files = []
    for file in files:
        if not os.path.exists(file):
            missing_files.append(file)
    
    if missing_files:
        print(f"缺失文件: {missing_files}")
        return False
    
    print("所有必需文件存在")
    print("\n安装完成!")
    print("\n使用方法:")
    print("python optimized_gad_mine.py -g protein_sequences.fasta -o results")
    print("\n更多信息请参考 README.md")
    
    return True

if __name__ == "__main__":
    install()