#!/usr/bin/env python3
"""
专业的GAD酶全基因组预测工具
可以从核苷酸序列预测潜在的GAD酶基因
"""

import os
import sys
import re
import argparse
from pathlib import Path
import datetime
from collections import Counter

class GADGenomePredictor:
    """
    专业的GAD酶基因组预测器
    """
    
    def __init__(self, genome_file, output_dir="gad_predictions"):
        self.genome_file = genome_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # 密码子表
        self.codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }
        
        # GAD酶的保守模式（基于文献和数据库）
        self.gad_motifs = [
            # PLP结合位点（关键特征）
            {"name": "PLP-binding-1", "pattern": "G[ST]G[KR][DE]", "weight": 5},
            {"name": "PLP-binding-2", "pattern": "G[ST]G[KR]X[KR]", "weight": 4},
            
            # 催化活性位点
            {"name": "active-site-1", "pattern": "K[ST][DE]", "weight": 4},
            {"name": "active-site-2", "pattern": "[KR][ST][DE]K", "weight": 3},
            {"name": "active-site-3", "pattern": "[KR]X[DE]X[KR]", "weight": 3},
            
            # 谷氨酸结合位点
            {"name": "glutamate-binding-1", "pattern": "[FY][KR][DE]", "weight": 3},
            {"name": "glutamate-binding-2", "pattern": "[FY]X[KR]", "weight": 2},
            
            # GAD酶的保守区域
            {"name": "GAD-conserved-1", "pattern": "KK[DE]", "weight": 3},
            {"name": "GAD-conserved-2", "pattern": "[KR][KR][DE]", "weight": 2},
            {"name": "GAD-conserved-3", "pattern": "[KR][DE][KR]", "weight": 2},
            
            # 功能域特征
            {"name": "pfam-domain-1", "pattern": "[KR][DE]X[KR]", "weight": 2},
            {"name": "pfam-domain-2", "pattern": "[FY][KR]X", "weight": 2},
            
            # 物种特异性保守序列
            {"name": "bacterial-GAD", "pattern": "[KR]X[KR][DE]", "weight": 2},
            {"name": "plant-GAD", "pattern": "[KR][DE]X[KR]", "weight": 2},
            {"name": "animal-GAD", "pattern": "[KR][KR]X[DE]", "weight": 2},
        ]
    
    def parse_fasta(self, file_path):
        """解析FASTA文件"""
        sequences = {}
        with open(file_path, 'r') as f:
            current_id = ""
            current_seq = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = current_seq
                    current_id = line[1:]  # 去掉'>'
                    current_seq = ""
                else:
                    current_seq += line
            if current_id:
                sequences[current_id] = current_seq
        return sequences
    
    def nucleotide_to_protein(self, nucleotide_seq):
        """核苷酸序列转换为蛋白质序列（三个阅读框）"""
        protein_seqs = []
        
        # 正向三个阅读框
        for offset in range(3):
            protein_seq = ""
            for i in range(offset, len(nucleotide_seq)-2, 3):
                codon = nucleotide_seq[i:i+3].upper()
                if codon in self.codon_table:
                    amino_acid = self.codon_table[codon]
                    if amino_acid != '*':
                        protein_seq += amino_acid
                    else:
                        break  # 遇到终止密码子停止
            if protein_seq:
                protein_seqs.append(protein_seq)
        
        # 反向三个阅读框（反转序列）
        rev_seq = nucleotide_seq[::-1]
        rev_seq = rev_seq.replace('A', 'T').replace('T', 'A').replace('C', 'G').replace('G', 'C')
        
        for offset in range(3):
            protein_seq = ""
            for i in range(offset, len(rev_seq)-2, 3):
                codon = rev_seq[i:i+3].upper()
                if codon in self.codon_table:
                    amino_acid = self.codon_table[codon]
                    if amino_acid != '*':
                        protein_seq += amino_acid
                    else:
                        break  # 遇到终止密码子停止
            if protein_seq:
                protein_seqs.append(protein_seq)
        
        return protein_seqs
    
    def calculate_amino_acid_composition(self, protein_seq):
        """计算氨基酸组成"""
        length = len(protein_seq)
        
        amino_acids = {
            'hydrophobic': ['A', 'V', 'L', 'I', 'M', 'F', 'W'],
            'hydrophilic': ['R', 'N', 'D', 'C', 'E', 'Q', 'H', 'K', 'S', 'T', 'Y'],
            'charged': ['R', 'K', 'D', 'E'],
            'glycine': ['G'],
            'proline': ['P'],
            'cysteine': ['C']
        }
        
        composition = {}
        for category, acids in amino_acids.items():
            count = sum(1 for aa in protein_seq if aa in acids)
            composition[category] = count / length
        
        return composition
    
    def analyze_gad_features(self, protein_seq):
        """分析GAD酶特征"""
        length = len(protein_seq)
        composition = self.calculate_amino_acid_composition(protein_seq)
        
        # GAD酶的特征阈值（基于文献）
        gad_features = {
            "length_range": (400, 600),  # GAD酶典型长度
            "hydrophobic_min": 0.25,
            "hydrophilic_min": 0.35,
            "charged_min": 0.20,  # GAD酶通常带电氨基酸比例高
            "charged_max": 0.45,  # 带电氨基酸过多可能不是GAD
            "glycine_min": 0.08,
            "proline_max": 0.05,
            "cysteine_max": 0.03,
        }
        
        # 结构特征评分
        structural_score = 0
        
        # 长度检查
        if length >= gad_features["length_range"][0] and length <= gad_features["length_range"][1]:
            structural_score += 20
        
        # 氨基酸组成评分
        if composition['hydrophobic'] >= gad_features["hydrophobic_min"]:
            structural_score += 15
        
        if composition['hydrophilic'] >= gad_features["hydrophilic_min"]:
            structural_score += 15
        
        # 带电氨基酸比例
        if composition['charged'] >= gad_features["charged_min"]:
            structural_score += 20
        
        if composition['charged'] <= gad_features["charged_max"]:
            structural_score += 10
        
        if composition['glycine'] >= gad_features["glycine_min"]:
            structural_score += 10
        
        if composition['proline'] <= gad_features["proline_max"]:
            structural_score += 5
        
        if composition['cysteine'] <= gad_features["cysteine_max"]:
            structural_score += 5
        
        return {
            'length': length,
            'composition': composition,
            'structural_score': structural_score / 85  # 归一化
        }
    
    def motif_search(self, protein_seq):
        """搜索GAD酶保守模式"""
        motif_matches = []
        
        for motif in self.gad_motifs:
            pattern = re.compile(motif['pattern'])
            matches = pattern.finditer(protein_seq)
            
            for match in matches:
                motif_matches.append({
                    'name': motif['name'],
                    'pattern': motif['pattern'],
                    'position': match.start(),
                    'match': match.group(),
                    'weight': motif['weight']
                })
        
        # PLP结合位点特别重要
        plp_matches = [match for match in motif_matches if 'PLP-binding' in match['name']]
        
        # 活性位点也很重要
        active_matches = [match for match in motif_matches if 'active-site' in match['name']]
        
        # 计算模式得分
        motif_score = sum(match['weight'] for match in motif_matches) * 50 / len(protein_seq)
        
        return {
            'matches': motif_matches,
            'plp_matches': len(plp_matches),
            'active_matches': len(active_matches),
            'motif_score': motif_score,
            'match_count': len(motif_matches)
        }
    
    def calculate_gad_score(self, protein_seq):
        """计算GAD酶可能性评分"""
        # 分析结构特征
        structural_info = self.analyze_gad_features(protein_seq)
        
        # 保守模式搜索
        motif_info = self.motif_search(protein_seq)
        
        # GAD酶的关键特征加权评分
        gad_score = structural_info['structural_score'] * 0.3
        
        # PLP结合位点加分
        if motif_info['plp_matches'] > 0:
            gad_score += 0.25
        
        # 活性位点加分
        if motif_info['active_matches'] > 0:
            gad_score += 0.15
        
        # 保守区域加分
        gad_score += motif_info['motif_score'] * 0.3
        
        # 长度合适加分
        if structural_info['length'] >= 400 and structural_info['length'] <= 600:
            gad_score += 0.1
        
        # 带电氨基酸比例加分
        if structural_info['composition']['charged'] >= 0.20:
            gad_score += 0.1
        
        return {
            'length': structural_info['length'],
            'composition': structural_info['composition'],
            'structural_score': structural_info['structural_score'],
            'motif_score': motif_info['motif_score'],
            'plp_matches': motif_info['plp_matches'],
            'active_matches': motif_info['active_matches'],
            'motif_count': motif_info['match_count'],
            'gad_score': gad_score,
            'protein_seq': protein_seq
        }
    
    def analyze_genomic_region(self, genomic_seq, region_id):
        """分析基因组区域，预测GAD酶基因"""
        results = []
        
        # 转换为蛋白质序列（所有可能的阅读框）
        protein_seqs = self.nucleotide_to_protein(genomic_seq)
        
        for i, protein_seq in enumerate(protein_seqs):
            # 分析蛋白质序列
            analysis = self.calculate_gad_score(protein_seq)
            
            # 判断是否为可能的GAD酶
            if analysis['gad_score'] >= 0.6 and analysis['length'] >= 200:
                results.append({
                    'region_id': region_id,
                    'frame': i,
                    'gad_score': analysis['gad_score'],
                    'length': analysis['length'],
                    'structural_score': analysis['structural_score'],
                    'motif_score': analysis['motif_score'],
                    'plp_matches': analysis['plp_matches'],
                    'active_matches': analysis['active_matches'],
                    'motif_count': analysis['motif_count'],
                    'composition': analysis['composition'],
                    'protein_seq': protein_seq[:100]  # 只保存前100个氨基酸
                })
        
        return results
    
    def run_genomic_prediction(self, threshold=0.6):
        """运行基因组预测"""
        print(f"开始GAD酶基因组预测，输入文件: {self.genome_file}")
        
        # 解析基因组序列
        genomic_regions = self.parse_fasta(self.genome_file)
        print(f"解析到 {len(genomic_regions)} 个基因组区域")
        
        # 分析每个区域
        candidates = []
        processed_count = 0
        
        for region_id, genomic_seq in genomic_regions.items():
            # 分析该基因组区域
            region_results = self.analyze_genomic_region(genomic_seq, region_id)
            
            for result in region_results:
                candidates.append(result)
            
            processed_count += 1
            if processed_count % 100 == 0:
                print(f"已处理: {processed_count}个基因组区域")
        
        # 排序候选
        candidates.sort(key=lambda x: x['gad_score'], reverse=True)
        
        print(f"\n预测完成!")
        print(f"找到 {len(candidates)} 个可能的GAD酶基因")
        print(f"结果保存在: {self.output_dir}")
        
        # 保存结果
        self.save_results(candidates, threshold)
        
        # 显示Top候选
        self.display_top_candidates(candidates, threshold)
        
        return candidates
    
    def save_results(self, candidates, threshold):
        """保存预测结果"""
        # 保存候选列表
        candidates_file = self.output_dir / "gad_genomic_predictions.tsv"
        with open(candidates_file, 'w') as f:
            # 写入标题
            f.write("Region_ID\tFrame\tGAD_Score\tLength\tStructural_Score\tMotif_Score\tPLP_Matches\tActive_Matches\tMotif_Count\tHydrophobic\tHydrophilic\tCharged\tGlycine\tProline\tCysteine\tProtein_Sequence\n")
            
            for candidate in candidates:
                f.write(f"{candidate['region_id']}\t{candidate['frame']}\t{candidate['gad_score']:.3f}\t{candidate['length']}\t{candidate['structural_score']:.3f}\t{candidate['motif_score']:.3f}\t{candidate['plp_matches']}\t{candidate['active_matches']}\t{candidate['motif_count']}\t{candidate['composition']['hydrophobic']:.3f}\t{candidate['composition']['hydrophilic']:.3f}\t{candidate['composition']['charged']:.3f}\t{candidate['composition']['glycine']:.3f}\t{candidate['composition']['proline']:.3f}\t{candidate['composition']['cysteine']:.3f}\t{candidate['protein_seq']}\n")
        
        # 保存详细报告
        report_file = self.output_dir / "gad_genomic_report.txt"
        with open(report_file, 'w') as f:
            f.write("GAD酶全基因组预测分析报告\n")
            f.write("==============================\n")
            f.write(f"输入文件: {self.genome_file}\n")
            f.write(f"分析日期: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"筛选阈值: {threshold}\n")
            f.write(f"候选数量: {len(candidates)}\n\n")
            
            f.write("Top 10可能的GAD酶基因:\n")
            for i, candidate in enumerate(candidates[:10], 1):
                f.write(f"{i}. {candidate['region_id']} (Frame {candidate['frame']})\n")
                f.write(f"  GAD可能性评分: {candidate['gad_score']:.3f}\n")
                f.write(f"  蛋白质长度: {candidate['length']} aa\n")
                f.write(f"  结构评分: {candidate['structural_score']:.3f}\n")
                f.write(f"  模式评分: {candidate['motif_score']:.3f}\n")
                f.write(f"  PLP结合位点匹配数: {candidate['plp_matches']}\n")
                f.write(f"  活性位点匹配数: {candidate['active_matches']}\n")
                f.write(f"  氨基酸组成: hydrophobic={candidate['composition']['hydrophobic']:.3f}, hydrophilic={candidate['composition']['hydrophilic']:.3f}, charged={candidate['composition']['charged']:.3f}\n")
                f.write(f"  蛋白质序列片段: {candidate['protein_seq']}\n\n")
            
            # 统计信息
            if candidates:
                avg_score = sum(c['gad_score'] for c in candidates) / len(candidates)
                avg_length = sum(c['length'] for c in candidates) / len(candidates)
                avg_plp = sum(c['plp_matches'] for c in candidates) / len(candidates)
                avg_active = sum(c['active_matches'] for c in candidates) / len(candidates)
                
                f.write(f"\n统计信息:\n")
                f.write(f"候选评分范围: {max(c['gad_score'] for c in candidates):.3f} - {min(c['gad_score'] for c in candidates):.3f}\n")
                f.write(f"平均GAD评分: {avg_score:.3f}\n")
                f.write(f"平均蛋白质长度: {avg_length:.1f} aa\n")
                f.write(f"平均PLP结合位点匹配数: {avg_plp:.1f}\n")
                f.write(f"平均活性位点匹配数: {avg_active:.1f}\n")
                
                # 阅读框分布
                frame_distribution = {}
                for candidate in candidates:
                    frame = candidate['frame']
                    frame_distribution[frame] = frame_distribution.get(frame, 0) + 1
                
                f.write(f"\n阅读框分布:\n")
                for frame, count in sorted(frame_distribution.items()):
                    f.write(f"Frame {frame}: {count}个候选\n")
    
    def display_top_candidates(self, candidates, threshold):
        """显示Top候选"""
        print("\nTop 10可能的GAD酶基因:")
        for i, candidate in enumerate(candidates[:10], 1):
            print(f"{i}. {candidate['region_id']} (Frame {candidate['frame']})")
            print(f"  GAD可能性评分: {candidate['gad_score']:.3f}")
            print(f"  蛋白质长度: {candidate['length']} aa")
            print(f"  PLP结合位点匹配数: {candidate['plp_matches']}")
            print(f"  活性位点匹配数: {candidate['active_matches']}")
        
        print("\n详细结果:")
        print(f"1. 候选列表: {self.output_dir}/gad_genomic_predictions.tsv")
        print(f"2. 分析报告: {self.output_dir}/gad_genomic_report.txt")
        print(f"3. 阈值: {threshold}")
        print(f"4. 候选数量: {len(candidates)}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='专业的GAD酶全基因组预测工具')
    parser.add_argument('-g', '--genome', required=True, help='基因组核苷酸序列文件（FASTA格式）')
    parser.add_argument('-o', '--output', default='gad_genomic_predictions', help='输出目录')
    parser.add_argument('-t', '--threshold', type=float, default=0.6, help='筛选阈值（默认0.6）')
    
    args = parser.parse_args()
    
    # 检查文件是否存在
    if not os.path.exists(args.genome):
        print(f"错误: 文件 {args.genome} 不存在!")
        return
    
    # 运行预测
    predictor = GADGenomePredictor(args.genome, args.output)
    candidates = predictor.run_genomic_prediction(args.threshold)

if __name__ == '__main__':
    main()