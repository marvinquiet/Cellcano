import os, sys, glob, shutil, re
import gzip
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

from typing import List

def _set_ArchR(genome: str, 
        seed: int = 2022, thread: int = 1) -> None:
    '''
    Basic settings for running ArchR
    ---
    Input: 
        - genome: mm9, mm10, hg19, hg38
        - threads: number of threads used to run ArchR 
            set to be 1 to avoid parallel related problem
    '''
    ArchR = importr('ArchR') ## import ArchR package in python
    ## import environment for different threads
    importr('parallel') 

    ArchR_thread_func = robjects.r['addArchRThreads']
    ArchR_thread_func(threads=thread)
    seed_func = robjects.r['set.seed']
    seed_func(seed)

    ## if the genome is not installed by BiocManager, ArchR will install it
    ArchR_genome_func = robjects.r['addArchRGenome']
    ArchR_genome_func(genome)

def _run_ArchR(input_files: List[str], sample_names: List[str], 
        output_dir: str = 'ArchROutput',
        min_tss: int = 4, min_frags: int = 1000,
        save_proj : bool = True):
    '''
    Run ArchR on input files
    ---
    Input: 
        - input_files: list of input files
        - sample_names: corresponding input sample names
        - min_tss: TSS enrichment score threshold to filter low quality cells
        - min_frags: number of mapped fragment threshold to filter low quality cells
        - save_proj: save ArchR project as RDS

    '''
    if len(input_files) != len(sample_names):
        ## exit the program
        sys.exit("Length of input files and sample names are not equal. Please make sure these two arguments are corresponding to each other.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
 
    ArchR_create_func = robjects.r['createArrowFiles']
    ArchR_arrowfiles = ArchR_create_func(inputFiles=robjects.StrVector(input_files), 
            sampleNames=robjects.StrVector(sample_names),
            minTSS=min_tss, minFrags=min_frags, 
            addGeneScoreMat=True)
    if len(ArchR_arrowfiles) != len(input_files):
        print("")
        ## TODO: add warning
        pass

    ArchR_proj_func = robjects.r['ArchRProject']
    ArchR_proj = ArchR_proj_func(ArrowFiles=robjects.StrVector(ArchR_arrowfiles),
            outputDirectory=output_dir, copyArrows=True) ## copy everything there and delete all files generated
    if save_proj:
        ArchR_save_func = robjects.r['saveArchRProject']
        ArchR_save_func(ArchRProj=ArchR_proj, outputDirectory=output_dir, load=False)

    ArchR_genescore_func = robjects.r['getMatrixFromProject']
    ArchR_genescore_mat = ArchR_genescore_func(ArchR_proj, useMatrix="GeneScoreMatrix")

    mm, base = importr('Matrix'), importr("base") ## import Matrix package in python
    dollar = base.__dict__["$"]
    mm_write_func, write_func = robjects.r['writeMM'], robjects.r['write']
    assay_func, rowData_func, colnames_func = robjects.r['assays'], robjects.r['rowData'], robjects.r['colnames']
    mm_write_func(dollar(assay_func(ArchR_genescore_mat), 'GeneScoreMatrix'),
            output_dir+os.sep+'ArchR_genescore.mtx')
    write_func(dollar(rowData_func(ArchR_genescore_mat), 'name'),
            output_dir+os.sep+'ArchR_genescore_genes.tsv')
    write_func(colnames_func(ArchR_genescore_mat),
            output_dir+os.sep+'ArchR_genescore_barcodes.tsv')

    genescore_path = output_dir+os.sep+'ArchR_genescore.mtx'
    if not os.path.exists(genescore_path):
        print("Error in producing gene score matrix!")
    else:
        ## gzip the gene score matrix
        with open(genescore_path, 'rb') as f_in:
            with gzip.open(genescore_path+'.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(genescore_path)

    ## deal with files
    cur_dir = os.getcwd()
    qc_dir = cur_dir+os.sep+'QualityControl'
    if os.path.exists(qc_dir):
        shutil.move(qc_dir, output_dir+os.sep+'QualityControl') ## check filter
    log_dir = cur_dir+os.sep+'ArchRLogs'
    if os.path.exists(log_dir):
        shutil.move(log_dir, output_dir+os.sep+'ArchRLogs') ## check log for ArchR
    rplot_file = cur_dir+os.sep+'Rplots.pdf'
    if os.path.exists(rplot_file):
        os.remove(rplot_file)

    tmp_dir = cur_dir+os.sep+'tmp'
    if os.path.exists(tmp_dir) and len(os.listdir(tmp_dir)) == 0:
        os.rmdir(tmp_dir)  ## remove tmp directories
    else:
        ## TODO: add warining
        pass

    arrowfiles_dir = output_dir+os.sep+'ArrowFiles'
    if os.path.exists(arrowfiles_dir) and len(os.listdir(arrowfiles_dir)) != 0:
        cur_arrowfiles = glob.glob(cur_dir+os.sep+'*.arrow')
        for x in cur_arrowfiles: os.remove(x)  ## remove arrow files
    else:
        ## TODO: add warning
        pass

def preprocess(args):
    '''Preprocess pipeline
    '''
    if not os.path.exists(args.input_dir):
        sys.exit("Error: %s does not exist." % args.input_dir)

    input_files = glob.glob(args.input_dir+os.sep+'*fragments.tsv.gz')
    input_bam_files = glob.glob(args.input_dir+os.sep+'*.bam')
    input_files.extend(input_bam_files)

    if len(input_files) == 0:
        sys.exit("Error: there are no files in the directory matching the patterns of: *fragments.tsv.gz, *.bam.")
    print("Following files will be processed: \n %s" % '\n\t'.join(input_files))

    if not args.sample_names:
        sample_names = []
        for input_file in input_files:
            if "fragments.tsv.gz" in input_file:
                sample_name = re.sub('_*fragments.tsv.gz$', '', os.path.basename(input_file))
            elif "bam" in input_file:
                sample_name = re.sub('.bam$', '', os.path.basename(input_file))
            sample_names.append(sample_name)
    else:
        sample_names = args.sample_names

    if sample_names and len(input_files) != len(sample_names):
        sys.exit("Error: the length of input files and sample names are unequal.")
    print("Following samples will be preprocessed: %s \n\n" % ' '.join(sample_names))

    if args.output_dir and os.path.exists(args.output_dir):
        output_dir = args.output_dir
    else:
        print("Output directory %s does not exist." % args.output_dir)
        print("Using input directory as output.\n\n")
        output_dir = args.input_dir
    print("Output will be written to %s" % output_dir)

    ## run ArchR pipeline
    _set_ArchR(genome=args.genome, thread=args.threads)
    _run_ArchR(input_files, sample_names, output_dir=output_dir)


if __name__ == '__main__':
    data_dir = "/home/wma36/gpu/Pyramid/Pyramid/samples"

    suffix = "_filtered_fragments.tsv.gz"
    input_files =  glob.glob(data_dir+os.sep+'*'+suffix)
    sample_names = ["PBMC_Rep1"]
    _set_ArchR(genome='hg19', seed=2022)
    _run_ArchR(input_files, sample_names, output_dir=data_dir, add_tilemat=True)

    #data_dir = "/home/wma36/gpu/data/Mousebrain_scATACseq/dscATACseq"

    #suffix = '.fragments.tsv.gz'
    #input_files =  glob.glob(data_dir+os.sep+'*'+suffix)
    #sample_names = [os.path.basename(x).replace(suffix, '') for x in input_files]

    #_set_ArchR(genome='mm10', seed=2022)
    #_run_ArchR(input_files, sample_names, output_dir=data_dir)

    #data_dir = "/home/wma36/gpu/data/humanPBMC_scATACseq/dscATACseq_isolated"
    #suffix = '.fragments.tsv.gz'
    #input_files = glob.glob(data_dir+os.sep+'*'+suffix)
    #sample_names = [os.path.basename(x).replace(suffix, '') for x in input_files]

    #_set_ArchR(genome='hg19', seed=2022)
    #_run_ArchR(input_files, sample_names, output_dir=data_dir, add_tilemat=True)

    #data_dir = "/home/wma36/gpu/data/humanPBMC_scATACseq/GSE139369_PBMC_MPAL_scATACseq"
    #suffix = '.fragments.tsv.gz'
    #input_files = glob.glob(data_dir+os.sep+'*'+suffix)
    #sample_names = [os.path.basename(x).replace(suffix, '') for x in input_files]

    #_set_ArchR(genome='hg19', seed=2022)
    #_run_ArchR(input_files, sample_names, output_dir=data_dir)

    #data_dir = "/home/wma36/gpu/data/Mousebrain_scATACseq/GSE111586_mouseatlas"
    #suffix = ".mm10.sorted.withbarcodes.bam"
    #input_files = glob.glob(data_dir+os.sep+'*'+suffix)
    #sample_names = [os.path.basename(x).replace(suffix, '') for x in input_files]

    #_set_ArchR(genome='mm10', seed=2022)
    #_run_ArchR(input_files, sample_names, output_dir=data_dir, min_tss=1, min_frags=100)

    ## human PBMC COVID
    #data_dir = "/home/wma36/gpu/data/humanPBMC_scATACseq/YangChen_COVID19"
    #suffix = "-fragments.tsv.gz"
    #input_files = glob.glob(data_dir+os.sep+'*'+suffix)
    #sample_names = [os.path.basename(x).replace(suffix, '').replace('500_', '') for x in input_files]

    #_set_ArchR(genome='hg19', seed=2022)
    #_run_ArchR(input_files, sample_names, output_dir=data_dir, min_tss=1, min_frags=1000)

