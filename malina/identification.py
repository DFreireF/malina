from ROOT import *
from pysimtof.importdata import *

class Identify():
    def __init__(self, data_dict): #data_dict= self.histogram_dict[:][1]
        histos=data_dict[:][0]
        data=data_dict[:][1]
        
    def tominimize(self, Brho):  
        self.calculate_ion_parameters(Brho)
        SRF = self.Frequence_rel*harm*self.SRRF[self.aux]
        tominimize = abs((self.ff[self.pp.argmax()]+self.fcenter)-SRF)
        print(f'tominimize={tominimize}, BRho={BRho}')
        return tominimize

    # Performs minimization of f_data[IsochroIon]-f_sample[RefIon(Brho)]
    def minimize_dist(self):
        print(f'function to minimize before minimizing: {self.tominimize(Brho)}')
        self.pdict['Brho'] = minimize(self.tominimize, [Brho], method='Powell',
                                      bounds=[(6.900, 6.910)], tol=1e-5).x[0]
        print(f'function minimized: {self.tominimize(Brho)}')
        self.calculate_ion_parameters(Brho)
        
        brho_correction = False
        if brho_correction == True:
            print(f"Brho initial: {self.pdict['Brho']}")
            self.BRhoCorrection()
            print(f"Brho final: {self.pdict['Brho']}")
            self.calculate_ion_parameters(Brho)
            self.SRF = [element*self.Frequence_Rel*Harmonic for element in self.SRRF]

    def set_peaks_exp(self):
        xpeaks=creategui.set_peaks('exp_data')
        #get_info_peaks(peaks(set_peaks))
        #get_pattern(this peaks with the others)--> distance matrix, etc...
        pattern=patternfinder.PatternFinder(data_dict['exp_data'][1][:,0], xpeaks)
        a, b =pattern.do_patter_matching()
        print(a,b)
        
if __name__='__main__':
    try:        
        mydata=ImportData(args.filename, args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge)
        mycanvas = CreateGUI(mydata.exp_data, mydata.simulated_data_dict, args.refisotope, mydata.nuclei_names)
        pattern=Identify(mycanvas.histogram_dict)
        pattern.set_peaks_exp()
        
    except: pass
