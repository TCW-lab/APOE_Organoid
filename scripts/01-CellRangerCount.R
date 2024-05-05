out<-'outputs/01-CellRangerCount'
dir.create(out)

source('../../utils/r_utils.R')

#based on : https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-gex-count
#synopsis
# cellranger count --id=sample345 \
#            --transcriptome=/opt/refdata-gex-GRCh38-2020-A \
#            --fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \ #path sep by coma if in different folder
#            --sample=mysample \ #name of the sample in the prefix of the fastq file
#            --localcores=8 \
#            --localmem=64

#test for 1
paste('cellranger count',
      '--id','S10',
      '--transcriptome','/projectnb/tcwlab/RefData/10XGenomics/refdata-gex-GRCh38-2020-A/',
      '--fastqs','/projectnb/tcwlab/RawData/RawRNAseq/scRNAseq/APOE_Jorganoid/data/data/230804_A00220_0594_BHG5H3DSX7/',
      '--sample','24')

#test for all

#get all sample id - sample name

samples_fastq<-list.files('/projectnb/tcwlab/RawData/RawRNAseq/scRNAseq/APOE_Jorganoid/data/data/230804_A00220_0594_BHG5H3DSX7/',pattern = '.fastq.gz$')
samples_dt<-unique(data.table(lib_id=str_extract(samples_fastq,'^[0-9]+'),
                       sample_id=str_extract(samples_fastq,'S[0-9]+')))

cmds<-sapply(samples_dt$sample_id, function(id){
  cmd<-paste('cellranger count',
             '--id',id,
             '--transcriptome','/projectnb/tcwlab/RefData/10XGenomics/refdata-gex-GRCh38-2020-A/',
             '--fastqs','/projectnb/tcwlab/RawData/RawRNAseq/scRNAseq/APOE_Jorganoid/data/data/230804_A00220_0594_BHG5H3DSX7/',
             '--sample',samples_dt[sample_id==id]$lib_id,
             '--output-dir',fp(out,id),
             '--localcores',16,
             '--localmem',128)
  return(cmd)
  
})
for(sample in names(cmds)){
  qsub_file<-fp('scripts',ps('01-CellRangerCount_',sample,'.qsub'))
  CreateJobFile(cmds[sample],
                file = qsub_file,
                nThreads = 16,maxHours = 48,
                modules = c('bcl2fastq/2.20','cellranger/7.2.0'))
  RunQsub(qsub_file,job_name = ps('CRcount',sample))
  
}


