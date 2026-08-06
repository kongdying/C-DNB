import sys
import os
target_path = os.path.abspath("../../my_libs")
if target_path not in sys.path:
    sys.path.insert(0, target_path)
import pandas as pd
import cmapPy
import gseapy as gp  
from cmapPy.pandasGEXpress.parse_gctx import parse
from concurrent.futures import ProcessPoolExecutor
df = pd.read_csv("./TF_T1.txt", sep="\t")
all_genes = df.iloc[:, 1].tolist()
up_list = [str(int(i)) for i in all_genes if pd.notna(i)]
down_list = []
geneset = {
    'up': up_list,
    'down': down_list
}
def GSEA_run_worker(batch_info):
    start, end = batch_info
    output_dir = "cmap/output_t1"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    outname = os.path.join(output_dir, f"result_{start}_{end}.txt")
    
    try:
        cidx = list(range(int(start), int(end)))
        gctx_obj = parse("../../data/level5_beta_trt_cp_n720216x12328.gctx", cidx=cidx)
        data_df = gctx_obj.data_df
        data_df.index = data_df.index.astype(str).str.replace(r'\.0$', '', regex=True)
    except Exception as e:
        print(f"Error reading GCTX {start}-{end}: {e}")
        return

    for sample_name in data_df.columns:
        try:
            # 排序：从高到低
            rnk = data_df[sample_name].sort_values(ascending=False)
            pre_res = gp.prerank(rnk=rnk,
                                 gene_sets=geneset,
                                 threads=1, # 进程内不再多线程防止资源竞争
                                 min_size=1,
                                 max_size=5000,
                                 permutation_num=100,
                                 outdir=None,
                                 seed=6,
                                 verbose=False)
            
            if not pre_res.res2d.empty:
                out = pre_res.res2d.copy()
                out.insert(0, 'cmap_sample_id', str(sample_name))
                out.to_csv(outname, sep="\t", index=False, header=False, mode="a")
        except Exception as e:
            continue
    print(f"Finished batch {start}-{end}")

if __name__ == '__main__':
    # 参数设置
    total_samples = 720216
    batch_size = 2000 # 并行时减小 batch_size 防止内存溢出
    max_workers = 8   # 根据你的服务器核心数调整
    
    # 构建批次列表
    batches = []
    for s in range(0, total_samples, batch_size):
        e = min(s + batch_size, total_samples)
        batches.append((s, e))

    print(f"开始并行计算，核心数: {max_workers}")
    
    # 调用进程池
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(GSEA_run_worker, batches)
df = pd.read_csv("./TF_T3.txt", sep="\t")
all_genes = df.iloc[:, 1].tolist()
up_list = [str(int(i)) for i in all_genes if pd.notna(i)]
down_list = []
geneset = {
    'up': up_list,
    'down': down_list
}
def GSEA_run_worker(batch_info):
    start, end = batch_info
    output_dir = "cmap/output_t3"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    outname = os.path.join(output_dir, f"result_{start}_{end}.txt")
    
    try:
        cidx = list(range(int(start), int(end)))
        gctx_obj = parse("../../data/level5_beta_trt_cp_n720216x12328.gctx", cidx=cidx)
        data_df = gctx_obj.data_df
        data_df.index = data_df.index.astype(str).str.replace(r'\.0$', '', regex=True)
    except Exception as e:
        print(f"Error reading GCTX {start}-{end}: {e}")
        return

    for sample_name in data_df.columns:
        try:
            # 排序：从高到低
            rnk = data_df[sample_name].sort_values(ascending=False)
            pre_res = gp.prerank(rnk=rnk,
                                 gene_sets=geneset,
                                 threads=1, # 进程内不再多线程防止资源竞争
                                 min_size=1,
                                 max_size=5000,
                                 permutation_num=100,
                                 outdir=None,
                                 seed=6,
                                 verbose=False)
            
            if not pre_res.res2d.empty:
                out = pre_res.res2d.copy()
                out.insert(0, 'cmap_sample_id', str(sample_name))
                out.to_csv(outname, sep="\t", index=False, header=False, mode="a")
        except Exception as e:
            continue
    print(f"Finished batch {start}-{end}")

if __name__ == '__main__':
    # 参数设置
    total_samples = 720216
    batch_size = 2000 # 并行时减小 batch_size 防止内存溢出
    max_workers = 8   # 根据你的服务器核心数调整
    
    # 构建批次列表
    batches = []
    for s in range(0, total_samples, batch_size):
        e = min(s + batch_size, total_samples)
        batches.append((s, e))

    print(f"开始并行计算，核心数: {max_workers}")
    
    # 调用进程池
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(GSEA_run_worker, batches)