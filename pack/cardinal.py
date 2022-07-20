from .core import Core
import os
import time


class Mapper(Core):
    def __init__(self, path, outputdir, reference, forks=4):
        super().__init__(path, outputdir, forks)
        self.reference = reference
        self.run()

    def run(self):
        for rgsm, fastgs in self.fastqs.items():
            R1, R2 = fastgs
            output = os.path.join(self.outputdir, f'{rgsm}.sam')
            out = os.path.join(self.outputdir, rgsm)
            getbwa = f'bwa mem -M -t {self.forks} -R "@RG\\tID:{rgsm}\\tSM:{rgsm}\\tLB:Genoks\\tPU:Illumina" {self.reference} {R1} {R2} > {output}'
            try:
                if os.system(getbwa) != 0:
                    raise Exception('Alignment rejected!')
            except Exception as e:
                print(e)
                print(
                    f'Couldnt mapp the {R1} and {R2} to the genome reference!')
                print('Going on in the line ...')
            else:
                print(f'Alignment for {R1} and {R2} finished successfully!')
                cmd = f'gatk --java-options "-Xms20g -Xmx30g" SortSam -I {output} -O {out}.sorted.bam -SO coordinate --CREATE_INDEX true'
                os.system(cmd)
                os.system(f'rm {output}')

class Sorter(Core):

    def __init__(self, path, outputdir, forks=4):
        super().__init__(path, outputdir, forks)
        self.func = self.getsort
        self.samples = self.sams
        self.multicore()

    def getsort(self, samplelist):
        rgsm, sam = samplelist
        out = os.path.join(self.outputdir, rgsm)
        #cmd = f'gatk --java-options "-Xms20g -Xmx30g" SortSam -I {sam} -O {out}.sorted.bam -SO coordinate --CREATE_INDEX true'
        #os.system(cmd)
        cmd = f'gatk --java-options "-Xms20g -Xmx30g" MarkDuplicates -I {out}.sorted.bam -O {out}.dedup.bam --METRICS_FILE {out}.markdup.metrics.txt --REMOVE_DUPLICATES true --CREATE_INDEX true'
        os.system(cmd)
        if os.path.isfile(f'{out}.dedup.bai'):
            cmd = f'rm {out}.sam {out}.sorted.bam {out}.sorted.bai'
            os.system(cmd)

# Read counts for


class ReadCount(Core):
    def __init__(self, path, outputdir, interval, reference, forks=4):
        super().__init__(path, outputdir, forks)
        self.interval = interval
        self.reference = reference
        self.func = self.function
        self.samples = self.bams
        self.multicore()

    def function(self, samplelist):
        rgsm, bam = samplelist
        out = os.path.join(self.outputdir, rgsm)
        cmd = f'gatk --java-options "-Xms20g -Xmx30g" CollectReadCounts -L {self.interval} -R {self.reference} -imr OVERLAPPING_ONLY -I {bam} -O {out}.TSV --format TSV'
        os.system(cmd)

# calculate coverage:
class Coverage(Core):

    def __init__(self, path, outputdir, reference, intervallist, bait_interval, forks=4):
        super().__init__(path, outputdir, forks)
        self.reference = reference
        self.bait_interval = bait_interval
        self.intervalList = intervallist
        self.func = self.function
        self.samples = self.bams
        self.multicore()

    def function(self, samplelist):
        rgsm, bam = samplelist
        out = os.path.join(self.outputdir, rgsm)
        cmd = f'gatk CollectHsMetrics -I {bam} -O {out}.coverage -R {self.reference} -TI {self.intervalList} -BI {self.bait_interval} --NEAR_DISTANCE 100' #--PER_BASE_COVERAGE {out}.per_base_cov --PER_TARGET_COVERAGE {out}.per_target_cov
        os.system(cmd)

#Depth of coverage
class Depth(Core):

    def __init__(self, path, outputdir, reference, intervallist, forks=4):
        super().__init__(path, outputdir, forks)
        self.reference = reference
        self.intervalList = intervallist
        self.func = self.function
        self.samples = self.bams
        self.multicore()

    def function(self, samplelist):
        rgsm, bam = samplelist
        out = os.path.join(self.outputdir, rgsm)
        cmd = f'gatk DepthOfCoverage -I {bam} -O {out}.depth -R {self.reference} -L {self.intervalList} -imr OVERLAPPING_ONLY --print-base-counts true'
        os.system(cmd)
        cmd = f' gatk CollectMultipleMetrics -I {bam} -O {out}'
        os.system(cmd)

class GvcfCaller(Core):
    def __init__(self, path, outputdir, reference, interval, padding, forks=4):
        super().__init__(path, outputdir, forks)
        self.reference = reference
        self.interval = interval
        self.func = self.function
        self.samples = self.bams
        self.padding = padding
        self.multicore()

    def function(self, samplelist):
        rgsm, bam = samplelist
        out = os.path.join(self.outputdir, rgsm)
        cmd = f'gatk --java-options "-Xms20g -Xmx30g" HaplotypeCaller -emit-ref-confidence GVCF -L {self.interval} --interval-padding {self.padding} -R {self.reference} -I {bam} -O {out}.g.vcf.gz'
        os.system(cmd)


class GetVCF(Core):
    def __init__(self, path, outputdir, reference, forks=4):
        super().__init__(path, outputdir, forks)
        self.reference = reference
        self.func = self.function
        self.samples = self.gvcfs
        self.multicore()

    def function(self, samplelist):
        rgsm, gvcf = samplelist
        out = os.path.join(self.outputdir, rgsm)
        cmd = f'gatk --java-options "-Xms20g -Xmx30g" GenotypeGVCFs -R {self.reference} -V {gvcf} -O {out}.vcf.gz'
        os.system(cmd)
        time.sleep(5)
        INdel = f'gatk --java-options "-Xms20g -Xmx30g" SelectVariants -R {self.reference} -V {out}.vcf.gz --select-type-to-include "INDEL" -O {out}.INDELs.vcf'
        Snp = f'gatk --java-options "-Xms20g -Xmx30g" SelectVariants -R {self.reference} -V {out}.vcf.gz --select-type-to-include "SNP" -O {out}.SNPs.vcf'
        os.system(INdel)
        os.system(Snp)
        time.sleep(5)
        IdelFilter = f'gatk --java-options "-Xms20g -Xmx30g" VariantFiltration -R {self.reference} -V {out}.INDELs.vcf -O {out}.FilIndel.vcf -filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "indel_filter" '
        SnpFilter = f'gatk --java-options "-Xms20g -Xmx30g" VariantFiltration -R {self.reference} -V {out}.SNPs.vcf -O {out}.FilSnp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "snp_filter"'
        os.system(IdelFilter)
        os.system(SnpFilter)
        time.sleep(5)
        combine = f'gatk --java-options "-Xms20g -Xmx30g" MergeVcfs -I {out}.FilIndel.vcf -I {out}.FilSnp.vcf -O {out}.hard_filtered.vcf.gz'
        os.system(combine)
        rm = f'rm {out}.FilIndel.vcf {out}.FilIndel.vcf.idx {out}.FilSnp.vcf {out}.FilSnp.vcf.idx {out}.INDELs.vcf {out}.INDELs.vcf.idx {out}.SNPs.vcf {out}.SNPs.vcf.idx'
        os.system(rm)

class VCFmetrics(Core):
    
    def __init__(self, path, outputdir, dbsnp, forks=4):
        super().__init__(path, outputdir, forks)
        self.func = self.function
        self.samples = self.gvcfs
        self.dbsnp = dbsnp
        self.multicore()
    
    def function(self, samplelist):
        rgsm, gvcf = samplelist
        out = os.path.join(self.outputdir, rgsm)
        cmd = f'gatk --java-options "-Xms20g -Xmx30g" CollectVariantCallingMetrics  -I {gvcf} -O {out}.VariantCallingMetrics --DBSNP {self.dbsnp} --GVCF_INPUT true'
        os.system(cmd)
