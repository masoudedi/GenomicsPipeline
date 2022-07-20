import os
import pandas as pd

class GatherMetrics:
    def __init__(self, path):
        self.path = path
        self.samples = []
        self.sampless = {}
        self.coverages = {}
        self.CovPropor = {}
        self.duplicates = {}
        self.lister()
        self.readCoverage()
        self.CovProportion()
        self.ReadDups()
        self.Mother()
        
    def lister(self):
        
        for i in os.listdir(self.path):
            if i.endswith('.coverage'):
                rgsm = os.path.basename(i).split('.')[0]
                self.samples.append(rgsm)
                self.sampless[rgsm] = {}
    
    def readCoverage(self):
        requested = ['ON_BAIT_BASES', 'NEAR_BAIT_BASES', 'OFF_BAIT_BASES', 'PCT_SELECTED_BASES', 'PCT_OFF_BAIT', 'ON_BAIT_VS_SELECTED', 'MEAN_BAIT_COVERAGE', 'PCT_USABLE_BASES_ON_BAIT', 'TOTAL_READS', 'PF_READS', 'PF_BASES', 'PF_BASES_ALIGNED', 'PF_UQ_BASES_ALIGNED', 'ON_TARGET_BASES', 'MEAN_TARGET_COVERAGE', 'MEDIAN_TARGET_COVERAGE', 'AT_DROPOUT', 'GC_DROPOUT']
        for sample in self.samples:
            count = 0
            fullpath = os.path.join(self.path, f'{sample}.coverage')
            with open(fullpath) as T:
                res = []
                for line in T:
                    if count == 6 or count == 7:
                        res.append(line.strip('\n').split('\t'))
                    count += 1
                res = {res[0][i]:res[1][i] for i in range(len(res[0]))}
                mmain = {}
                for key, value in res.items():
                    if key in requested:
                        mmain[key] = value
                # res = dict(list(res.items())[3:65])
                self.sampless[sample]['coverage'] = mmain
                self.coverages[sample] = mmain
    
    def CovProportion(self):
        for sample in self.samples:
            fullpath = os.path.join(self.path, f'{sample}.depth.sample_cumulative_coverage_proportions')
            with open(fullpath) as T:       
                res = []
                for line in T:
                    res.append(line.strip('\n').split(','))
                res = {res[0][i]:res[1][i] for i in range(len(res[0]))}
                res = dict(list(res.items())[1:102])
                self.sampless[sample]['CovPropor'] = res
                self.CovPropor[sample] = res

    def ReadDups(self):
        for sample in self.samples:
            count = 0
            fullpath = os.path.join(self.path, f'{sample}.markdup.metrics.txt')
            with open(fullpath) as T:
                res = []
                for line in T:
                    if count == 6 or count == 7:
                        res.append(line.strip('\n').split('\t'))
                    count += 1
                res = {res[0][i]:res[1][i] for i in range(len(res[0]))}
                self.sampless[sample]['duplicates'] = res
                self.duplicates[sample] = res
    
    def Mother(self):
        result = []
        #get dinamic header
        if self.samples:
            tt = self.sampless[self.samples[0]]
            header = ['Sample_name'] + list(tt['coverage'].keys()) + list(tt['duplicates'].keys()) + list(tt['CovPropor'].keys())
            result.append(header)
        for sample in self.sampless:
            tt = self.sampless[sample]
            resss = [sample] + list(tt['coverage'].values()) + list(tt['duplicates'].values()) + list(tt['CovPropor'].values())
            result.append(resss)
        df = pd.DataFrame(result)
        with pd.ExcelWriter(os.path.join(self.path, 'Metrics.Ilyome.xlsx')) as writer:
            df.to_excel(writer, sheet_name='Metrics', index=False, header=False)
