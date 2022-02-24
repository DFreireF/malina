from ROOT import *
from iqtools import *

class Identify_Isomer():
    def __init__(self, filename): #data_dict= self.histogram_dict[:][1]
        
    @staticmethod    
    def _cut_spectrogram(spec_to_cut, f1, f2, y1=None, y2=None):
        xcen=(f1+f2-2*self.iq.center)/2
        xspan=abs(f1-f2)
        nxx, nzz, nyy= get_cut_spectrogram(spect_to_cut[0], spec_to_cut[1], spec_to_cut[2], xcen=xcen, \
xspan=xspan)
        return (np.stack((nxx, nzz, nyy), axis=1)).reshape((len(nxx),3))

    def _create_spectrogram(self, howlong_t, where_t, lf, method=None):
        nframes=int(howlong_t*self.iq.fs/lf)
        sframes=int(where_t*self.iq.fs/lf)
        self.iq.read(nframes=nframes, lframes=lf, sframes=sframes)
        iq.method='mtm' #'fft', 'mtm', 'welch'                                                          
        if method: iq.method = method
        xx, yy, zz = self.iq.get_spectrogram(nframes,lframes)
        return (np.stack((xx, zz, yy), axis=1)).reshape((len(xx),3))
    
    @staticmethod
    def method1(filename, lframes, inject_t):
        self.iq = get_iq_object(filename)
        self.iq.read_samples(1)
        time=inject_t-2
        skip_time=0
        bg1 = self._create_spectrogram(time, skip_time, lframes)
        nbg1= self._cut_spectrogram(bg1, f1, f2)

        time=0.100 #100ms of data                                                                       
        skip_time=inject_t
        interest_data = self._create_spectrogram(time, skip_time, lframes)
        idata= self._cut_spectrogram(interest_data, f1, f2)

        time=inject_t-2
        skip_time=inject_t+1
        bg2 = ImportData._create_spectrogram(time, skip_time, lframes)
        nbg2= self._cut_spectrogram(bg2, f1, f2)

        self.background=


###########################
##########################
        
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
        
if __name__=='__main__':
    try:        
        mydata=ImportData(args.filename, args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge)
        mycanvas = CreateGUI(mydata.exp_data, mydata.simulated_data_dict, args.refisotope, mydata.nuclei_names)
        pattern=Identify(mycanvas.histogram_dict)
        pattern.set_peaks_exp()
        
    except: pass
