import os
import concurrent.futures as cf

class Core:
    def __init__(self, path, outputdir, forks=4):
        self.path = path
        self.outputdir = outputdir
        self.forks = forks
        self.func = None
        self.fastqs = {}
        self.bams = []
        self.sams = []
        self.gvcfs = []
        self.getlist()
        self.samples = None

    def getlist(self):
        for i in os.listdir(self.path):
            #gets list of fastqs:
            if i.endswith('_R1_001.fastq.gz'):
                rgsm = os.path.basename(i)[:-16]
                R2 = os.path.join(self.path, f'{rgsm}_R2_001.fastq.gz')
                R1 = os.path.join(self.path, f'{rgsm}_R1_001.fastq.gz')
                if os.path.isfile(R2):
                    self.fastqs[rgsm] = [R1, R2]
            
            if i.endswith('_1.fq.gz'):
                rgsm = os.path.basename(i)[:-8]
                R2 = os.path.join(self.path, f'{rgsm}_2.fq.gz')
                R1 = os.path.join(self.path, f'{rgsm}_1.fq.gz')
                if os.path.isfile(R2):
                    self.fastqs[rgsm] = [R1, R2]
            
            #starts getting list of sams:
            if i.endswith('sorted.bam'):
                rgsm = os.path.basename(i).split('.')[0]
                sam = os.path.join(self.path, os.path.basename(i))
                if os.path.isfile(sam):
                    self.sams.append([rgsm, sam])
            
            #starts getting list of bams
            if i.endswith('bam'):
                rgsm = os.path.basename(i).split('.')[0]
                bam = os.path.join(self.path, os.path.basename(i))
                if os.path.isfile(bam):
                    self.bams.append([rgsm, bam])
            
            #get list of g.vcfs
            if i.endswith('.g.vcf.gz'):
                rgsm = os.path.basename(i).split('.')[0]
                gvcf = os.path.join(self.path, os.path.basename(i))
                if os.path.isfile(gvcf):
                    self.gvcfs.append([rgsm, gvcf])
    
    def multicore(self):
        with cf.ProcessPoolExecutor(max_workers=40) as executor:
            workers = [executor.submit(self.func, samplelist)
                       for samplelist in self.samples]
            for p in cf.as_completed(workers):
	            p.result()
