#!/usr/bin/env python3
"""
优化版GAD酶基因组挖掘工具
"""

import os
import sys
import re
import json
import argparse
from pathlib import Path

class GADFinder:
    """
    GAD酶挖掘器
    """
    
    def __init__(self, genome_file, output_dir="results"):
        self.genome_file = genome_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # GAD酶的保守特征模式
        self.gad_motifs = [
            # PLP结合位点相关模式
            {"name": "PLP-binding", "pattern": "G[ST]G[KR][DE]", "weight": 3},
            {"name": "active-site", "pattern": "[KR][ST][DE]X[KR]", "weight": 2},
            {"name": "glutamate-binding", "pattern": "[FY]X[KR][DE]", "weight": 2},
            
            # GAD酶的典型氨基酸序列特征
            {"name": "conserved-GAD", "pattern": "KK[DE]K", "weight": 1},
            {"name": "catalytic-core", "pattern": "XX[KR]X[DE]", "weight": 1},
            
            # 基于已知GAD酶的保守区域
            {"name": "GAD-conserved-1", "pattern": "K[KR]X[DE]K", "weight": 1},
            {"name": "GAD-conserved-2", "pattern": "F[KR]X[DE]", "weight": 1},
        ]
        
        # GAD酶的参考序列（简化的示例序列）
        self.reference_gad_sequences = {
            "Ecoli_GadA": "MKYEPLIKKSPQKHFVLEGDNKNGEGYKQHSKWEGVEMYLKKVLMSCEAQPVNHSCVDLIKEVGEVAATGKKTGEGMVKKIRESIFGKPVPLTAGCGIMVSDNVKSIESLKLHKNLCPSIVVGGKKAKKDYPDCLKQALAVYATNFEGICPTMTEQDEENNTYYWRGHDDQHFHNEKEKLEDLTWTVKQRTQQKEITVRFYTNMNITIMAKVPYENFTSAVKEINEENVQFYKKDLEAVKRHQFVDCYQLYMSTQLQSAFGDKFPTFGVIGGISASVRHCLKECIKKIMEKDEAQFHNHQMNCSMPQCQGIGNVSLITWIPQRSDPLVTRAPKKDEKYLCETITPFKPDKLACLADMGPHYTLEVTKNGFNQIHMKLPGKNETPYLKIDNVQTKQTVGKNFNICLLVMAKQSNYQGIKVQYSNFICSNMDEAFGKMIKYNGKPFFSFKPVSYKPIKYTNIIGILKQVKEIDRFKDEKIEQRIVEKVEKKQNHFNYDQIEHFRQQHVGEYKNIYEKAKDMYEVYNHTGDYPVKGINPESQSCVDFEHAIQVCHSEHLCDHFNKHIFVVHRPASEKKAFTSKYYKKGYDTIENQIVVKLNEFETEYNIKHGFSHIQRQDVEKDYYREGNWDYHIFIKHEVDKECQYVLNCVTRNKVGIKNDNNQYSYINKIYSKDFSHFKQVKYLK",
            "Lactobacillus_Gad": "MKIKKIVIVKGRKSGEGVEEKNKKASFVVKAIVPNGTEIVQVKKRGSCVTFKDTYHGGFSVETNQRIKKSLNYTYFNGVVNDIKPVEGKYNKFVQDQYAELTKADYYVMCKLGVELFSLKKIDVGYEMNKKQYGEVCASLCETKQKLIYMVIPKDIYKDDHDFGDYEVKPDQVCRAAHQKMDIFIQGNEAYQDLTKNYYVHYKGQYSVWRNNKELEQVSCQLSLGDYDNKVQKFFENLNQQAHNFLYDKFIEIRNLDYQKFGKKKKYQKHRSVYKGGSLQDCQVTYCEQLLNRKQY",
            "Human_GAD67": "MMIQQIRANQAVQGFSLDLQVFDDCTKYATSMVEKGLMIKLNTYFNLFRAKPMSLIESPRIEAVTIVRFDNQCNFYQDSLWKWSCDSRYSLEGFDIVQIEIQGARDPQFCLSVAFKELGTTHPNIRUQLVIVKDHYSKHDHTSQFEDIYRRAFLRKWPYKFEPCTNIDEYTTAMVTMSKGPKSLYASJYTAPGWASFDVSQPECYLQYVPAIRKLKEHQLMNYSRTREPILQVSCYLMRKEQTYLKFSIPSHGFQPRIVWVLSADTYIYJYLPOKCYTDXFVJPDHNSHAYEMTQDLKPFLKJHASDTIRYWLIEKHERDNQSHFKHLKTRVTYNKCMLNRQNNSNLKELJHYGGTSFNHQGRNDIVFTQHEDSTQFFHYEIHNRQIRLQGSNSKGTRYMINSDLHSQEYLFYEVTNLEPMLVVDLFEVFFEIQSYVNTMFKLNCIHNMLQUYFQFPRKMLSVREGMKYKRHDRKCHCHLMSEWMYTLJTLPTTTRVQLFHDCNDPDPSLQFDNFMKQAJMNRTQYTLJSNGTYPMQYMHJVYSJRHINQJNSLHEYSMJNEVJTGTMNLYRQRMKEHQRLTCGDJFTSMTETESQNNYKDYKVNKIMKFVLIGVILKDDMKTDYREFLMFHYSGYFJPGQVGGDPKGIVLSYFLTDQWPRRSFGHVEEVRQPNJMQSXLTQMKSMQSNKPVNHGGGTQQLGHMSSJHVVYSKLGDYDM",
        }
    
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
    
    def analyze_sequence(self, sequence):
        """分析序列特征"""
        length = len(sequence)
        
        # 氨基酸组成分析
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
            count = sum(1 for aa in sequence if aa in acids)
            composition[category] = count / length
        
        # GAD酶的特征阈值
        gad_features = {
            "length_range": (400, 650),
            "hydrophobic_min": 0.25,
            "hydrophilic_min": 0.35,
            "charged_min": 0.15,
            "glycine_min": 0.08,
            "proline_max": 0.05,
        }
        
        # 结构特征评分
        structural_score = 0
        
        # 长度检查
        if length >= gad_features["length_range"][0] and length <= gad_features["length_range"][1]:
            structural_score += 20
        
        # 氨基酸组成评分
        if composition['hydrophobic'] >= gad_features["hydrophobic_min"]:
            structural_score += 10
        
        if composition['hydrophilic'] >= gad_features["hydrophilic_min"]:
            structural_score += 10
        
        if composition['charged'] >= gad_features["charged_min"]:
            structural_score += 10
        
        if composition['glycine'] >= gad_features["glycine_min"]:
            structural_score += 10
        
        if composition['proline'] <= gad_features["proline_max"]:
            structural_score += 5
        
        # 归一化得分
        normalized_score = structural_score / 65
        
        return {
            'length': length,
            'composition': composition,
            'structural_score': normalized_score,
            'structural_raw_score': structural_score
        }
    
    def motif_search(self, sequence):
        """搜索GAD酶保守模式"""
        motif_matches = []
        
        for motif in self.gad_motifs:
            pattern = re.compile(motif['pattern'])
            matches = pattern.finditer(sequence)
            
            for match in matches:
                motif_matches.append({
                    'name': motif['name'],
                    'pattern': motif['pattern'],
                    'position': match.start(),
                    'match': match.group(),
                    'weight': motif['weight']
                })
        
        # 计算模式匹配得分
        motif_score = sum(match['weight'] for match in motif_matches) * 100 / len(sequence)
        
        return {
            'matches': motif_matches,
            'motif_score': motif_score,
            'match_count': len(motif_matches)
        }
    
    def homology_comparison(self, sequence):
        """与参考序列进行同源性比较"""
        similarity_scores = []
        
        for ref_name, ref_seq in self.reference_gad_sequences.items():
            # 计算简单的相似性分数
            # 基于共享模式的数量
            seq_motifs = self.motif_search(sequence)
            ref_motifs = self.motif_search(ref_seq)
            
            shared_motifs = 0
            for seq_match in seq_motifs['matches']:
                for ref_match in ref_motifs['matches']:
                    if seq_match['name'] == ref_match['name']:
                        # 考虑模式位置的相似性
                        position_diff = abs(seq_match['position'] - ref_match['position'])
                        if position_diff < 100:  # 允许一定位置偏移
                            shared_motifs += seq_match['weight']
            
            # 序列长度相似性
            len_diff = abs(len(sequence) - len(ref_seq))
            len_similarity = 1 - (len_diff / max(len(sequence), len(ref_seq)))
            
            # 氨基酸组成相似性
            seq_features = self.analyze_sequence(sequence)
            ref_features = self.analyze_sequence(ref_seq)
            
            comp_diff = sum([
                abs(seq_features['composition']['hydrophobic'] - ref_features['composition']['hydrophobic']),
                abs(seq_features['composition']['hydrophilic'] - ref_features['composition']['hydrophilic']),
                abs(seq_features['composition']['charged'] - ref_features['composition']['charged']),
                abs(seq_features['composition']['glycine'] - ref_features['composition']['glycine'])
            ])
            
            composition_similarity = 1 - comp_diff
            
            # 综合相似性评分
            similarity_score = (shared_motifs * 0.4 + len_similarity * 0.3 + composition_similarity * 0.3)
            
            similarity_scores.append({
                'ref_name': ref_name,
                'similarity_score': similarity_score,
                'shared_motifs': shared_motifs,
                'len_similarity': len_similarity,
                'composition_similarity': composition_similarity
            })
        
        # 找到最佳匹配
        best_match = max(similarity_scores, key=lambda x: x['similarity_score'])
        
        return {
            'similarity_scores': similarity_scores,
            'best_match': best_match['ref_name'],
            'best_score': best_match['similarity_score']
        }
    
    def find_gad_candidates(self, protein_file):
        """在蛋白质文件中寻找GAD酶候选"""
        sequences = self.parse_fasta(protein_file)
        
        candidates = []
        
        for protein_id, sequence in sequences.items():
            print(f"分析: {protein_id}")
            
            # 1. 结构特征分析
            structural_features = self.analyze_sequence(sequence)
            
            # 2. 保守模式搜索
            motif_results = self.motif_search(sequence)
            
            # 3. 同源性比较
            homology_results = self.homology_comparison(sequence)
            
            # 综合评分
            total_score = (
                structural_features['structural_score'] * 0.4 +
                motif_results['motif_score'] * 0.3 +
                homology_results['best_score'] * 0.3
            )
            
            # 筛选阈值
            if total_score > 0.35:
                candidates.append({
                    'protein_id': protein_id,
                    'sequence_length': len(sequence),
                    'total_score': total_score,
                    'structural_score': structural_features['structural_score'],
                    'motif_score': motif_results['motif_score'],
                    'homology_score': homology_results['best_score'],
                    'best_match': homology_results['best_match'],
                    'motif_count': motif_results['match_count'],
                    'composition': structural_features['composition']
                })
        
        # 排序候选
        candidates.sort(key=lambda x: x['total_score'], reverse=True)
        
        return candidates
    
    def generate_report(self, candidates, output_file):
        """生成分析报告"""
        import time
        
        with open(output_file, 'w') as f:
            f.write("GAD酶挖掘分析报告\n")
            f.write("=======================\n")
            f.write(f"输入文件: {self.genome_file}\n")
            f.write(f"分析日期: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"候选数量: {len(candidates)}\n")
            f.write("\n")
            
            if candidates:
                f.write("Top 10候选GAD酶:\n")
                for i, candidate in enumerate(candidates[:10]):
                    f.write(f"{i+1}. {candidate['protein_id']}\n")
                    f.write(f"  综合评分: {candidate['total_score']:.3f}\n")
                    f.write(f"  结构评分: {candidate['structural_score']:.3f}\n")
                    f.write(f"  模式评分: {candidate['motif_score']:.3f}\n")
                    f.write(f"  同源性评分: {candidate['homology_score']:.3f}\n")
                    f.write(f"  最佳参考匹配: {candidate['best_match']}\n")
                    f.write(f"  匹配模式数量: {candidate['motif_count']}\n")
                    f.write(f"  氨基酸组成: hydrophobic={candidate['composition']['hydrophobic']:.3f}, hydrophilic={candidate['composition']['hydrophilic']:.3f}, charged={candidate['composition']['charged']:.3f}, glycine={candidate['composition']['glycine']:.3f}\n")
                    f.write("\n")
            else:
                f.write("未找到符合条件的候选GAD酶\n")
            
            # 统计信息
            f.write("\n统计信息:\n")
            if candidates:
                f.write(f"候选评分范围: {candidates[0]['total_score']:.3f} - {candidates[-1]['total_score']:.3f}\n")
                f.write(f"平均长度: {sum(c['sequence_length'] for c in candidates)/len(candidates)} aa\n")
                
                # 参考匹配分布
                match_distribution = {}
                for candidate in candidates:
                    match = candidate['best_match']
                    match_distribution[match] = match_distribution.get(match, 0) + 1
                
                f.write("\n参考序列匹配分布:\n")
                for match, count in match_distribution.items():
                    f.write(f"{match}: {count}个候选\n")
            else:
                f.write("没有候选序列\n")
    
    def run(self):
        """运行完整的GAD酶挖掘流程"""
        print(f"开始GAD酶挖掘，输入文件: {self.genome_file}")
        
        # 检查文件是否存在
        if not os.path.exists(self.genome_file):
            print(f"文件不存在: {self.genome_file}")
            return
        
        # 假设输入已经是蛋白质序列文件
        protein_file = self.genome_file
        
        # 寻找GAD酶候选
        candidates = self.find_gad_candidates(protein_file)
        
        # 保存结果
        results_file = self.output_dir / "gad_candidates.tsv"
        with open(results_file, 'w') as f:
            f.write("Protein_ID\tLength\tTotal_Score\tStructural_Score\tMotif_Score\tHomology_Score\tBest_Match\tMotif_Count\tHydrophobic\tHydrophilic\tCharged\tGlycine\tProline\tCysteine\n")
            
            for candidate in candidates:
                f.write(f"{candidate['protein_id']}\t{candidate['sequence_length']}\t"
                        f"{candidate['total_score']:.3f}\t{candidate['structural_score']:.3f}\t"
                        f"{candidate['motif_score']:.3f}\t{candidate['homology_score']:.3f}\t"
                        f"{candidate['best_match']}\t{candidate['motif_count']}\t"
                        f"{candidate['composition']['hydrophobic']:.3f}\t{candidate['composition']['hydrophilic']:.3f}\t"
                        f"{candidate['composition']['charged']:.3f}\t{candidate['composition']['glycine']:.3f}\t"
                        f"{candidate['composition']['proline']:.3f}\t{candidate['composition']['cysteine']:.3f}\n")
        
        # 生成报告
        report_file = self.output_dir / "gad_report.txt"
        self.generate_report(candidates, report_file)
        
        print(f"\n分析完成!")
        print(f"找到 {len(candidates)} 个候选GAD酶")
        print(f"结果保存在: {self.output_dir}")
        
        if candidates:
            print(f"\nTop 5候选:")
            for i, candidate in enumerate(candidates[:5]):
                print(f"{i+1}. {candidate['protein_id']}")
                print(f"  综合评分: {candidate['total_score']:.3f}")
                print(f"  最佳参考匹配: {candidate['best_match']}")
        
        return candidates


def main():
    import time
    
    parser = argparse.ArgumentParser(description='优化版GAD酶挖掘工具')
    parser.add_argument('-g', '--genome', required=True, help='蛋白质序列文件（FASTA格式）')
    parser.add_argument('-o', '--output', default='results', help='输出目录')
    parser.add_argument('-t', '--threshold', type=float, default=0.35, help='筛选阈值（默认0.35）')
    
    args = parser.parse_args()
    
    # 创建挖掘器
    miner = GADFinder(args.genome, args.output)
    
    # 运行挖掘
    candidates = miner.run()
    
    if candidates:
        print(f"\n详细结果:")
        print(f"1. 候选列表: {args.output}/gad_candidates.tsv")
        print(f"2. 分析报告: {args.output}/gad_report.txt")
        print(f"3. 阈值: {args.threshold}")
        print(f"4. 候选数量: {len(candidates)}")
    
    return 0


if __name__ == "__main__":
    main()