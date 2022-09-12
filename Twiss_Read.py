#This file reads PDF tables from the requested pages
import pymadx
import numpy as np
#This class reads in twiss files and converts them into G4BL gradients
class Twiss_read():
    def  __init__(self,file,ebtoggle=1):
        self.tw = pymadx.Data.Tfs(file)
        #print(self.tw)
        self.tw.ReportPopulations()
        self.tw.Plot()
        self.ebtoggle=ebtoggle #0 for B values, 1 for E values

        global m,E,k,p,l
        self.m = 938.2720813 # MeV c-2 mass
        self.Ek = 0.1 # MeV energy
        self.p = np.sqrt( (self.m+self.Ek)**2 - self.m**2 ) # MeV/c momentum
        self.p = self.p/1000 #convert to GeV/c
        self.l = 0.1
        self.element_size=10

    def k2g(self,k,l): #Converts from focusing strength into gradient, returns both a magnetic and electric strength.
        q=-1 #Charge of particle
        psi = self.p/0.299792458
        brel = self.p/(self.Ek+self.m)
        brel = 1.45987e-05
        # magnetic gradient
        g_m = q*(k*self.p)/(0.299792458) #T/m
        v = brel*2.99792458e8
        g_e = q*(k*(psi)*v)/1000000 # works /1000000 to convert V/m -> MV/m
        return g_m*l, g_e*l


    def getquads(self):	#Gets all quadrupoles based on the TFS keyword "QUADRUPOLE", adjust as required
        quads = []
        for el in self.tw:
            if el['KEYWORD'] == 'QUADRUPOLE' :
                k1l = el['K1L']
                ##print(el['NAME'],k1l)
                grad = self.k2g(k1l,self.element_size)
                grad = grad[self.ebtoggle]
                #cd = {el['NAME'] : "Name", grad:"grad",el['K1L']:"k1L"}
                s=float(el['S']) *1000
                cd = {"name": el['NAME'], "grad" : grad , "k1L" : el['K1L'],"S":s,"type" : el['KEYWORD']}
                quads.append(cd)
        return quads

    def getbends(self):	#Grabs all bends from TFS file with keyword "MATRIX", adjust as required!
        bends = []
        for el in self.tw:
            if el['KEYWORD'] == 'MATRIX' :	
                k1l = el['K1L']
                grad = self.k2g(k1l,self.element_size)
                grad = grad[0]
                s=float(el['S']) *1000
                cd = {"name": el['NAME'],"grad" : grad , "k1L" : el['K1L'],"S":s,"type" : el['KEYWORD']}
                bends.append(cd)
        return bends

    def getVkickers(self):	#Vertical Kickers
        Vkickers = []
        for el in self.tw:
            if el['KEYWORD'] == 'VKICKER' :
                k1l = el['K1L']
                grad = self.k2g(k1l,self.element_size)
                grad = grad[0]
                s=float(el['S']) *1000
                strname = el['NAME']
                strname.replace(".","_")
                cd = {"name": strname, "grad" : grad , "k1L" : el['K1L'],"S":s,"type" : el['KEYWORD']}
                Vkickers.append(cd)
        return Vkickers

    def getHkickers(self):	#Horizontal Kickers
        Hkickers=[]
        for el in self.tw:
            if el['KEYWORD'] == 'HKICKER' :
                k1l = el['K1L']
                grad = self.k2g(k1l,self.element_size)
                grad = grad[0]
                s=float(el['S']) *1000
                strname = el['NAME']
                strname.replace(".","_")
                cd = {"name": strname, "grad" : grad , "k1L" : el['K1L'],"S":s,"type" : el['KEYWORD']}
                Hkickers.append(cd)
        return Hkickers

    def getall(self):	#Gets all objects!
        all=[]
        for el in self.tw:
            if el['KEYWORD'] != 'DRIFT' and el['KEYWORD'] != 'MARKER':
                k1l = el['K1L']
                grad = self.k2g(k1l, self.element_size)
                #print(el['NAME'])
                if el['KEYWORD'] == 'QUADRUPOLE':
                    grad = grad[self.ebtoggle]
                else:
                    grad = grad[0]
                s=float(el['S']) *1000
                strname = el['NAME']
                strname=strname.replace(".","_")
                #print(strname)
                #cd = {"name": strname, "grad" : grad , "k1L" : el['K1L'],"S":s,"type" : el['KEYWORD']}
                cd = {"name": el['NAME'], "grad": grad, "k1L": el['K1L'], "S": s, "type": el['KEYWORD']}
                all.append(cd)
            elif (el['KEYWORD'] == 'MARKER'):
                     if el['NAME'] == "LNE04$END" or el['NAME'] == 'LNE02$END':
                         print("Using end point : %s \n" % (el['NAME']))
                         #print("FOUND END DETAILS: \n Location : %s" , (el['S']))
                         cd = {"name":el['NAME'],"S":el['S']*1000,"type":"EndLine04"}
                         all.append(cd)
        return all

    def getMakrers(self):
        markers=[]
        for el in self.tw:
            if el['KEYWORD'] == 'MARKER':
                markers.append(el)
        return markers

    def getBeamstart(self,beamname="LNE00LNE01LNE03LNE04$START"):	#Beam start - adjust either start beamname here or in whichever program you use to call these functions
        for el in self.tw:
            if el['NAME'] == beamname:
                sigmax = el['SIGMAX']
                sigmay = el['SIGMAY']
                dxbeta = el['DXBETA']
                dybeta = el['DYBETA']
                dpxbeta = el['DPXBETA']
                dpybeta = el['DPYBETA']
                sigmaxp = el['SIGMAXP']
                sigmayp = el['SIGMAYP']
                betx = el['BETX']
                S = el['S']
                X = el['X']
                Y = el['Y']
                DX = el['DX']
                PX = el['PX']
                DPX = el['DPX']
                DY = el['DY']
                PY = el['PY']
                DPY = el['DPY']
                #print
                tbeam = "beam gaussian nEvents=300 particle=anti_proton x=%s y=%s z=-200 sigmaX=%s sigmaXp=%s sigmaY=%s sigmaYp=%s" %(X,Y,sigmax/1000,sigmaxp,sigmay/1000,sigmayp)
                print(tbeam)
                return(el)

    def getkickersall(self):
        kickers=[]
        for el in self.tw:
            if el['KEYWORD'] == 'VKICKER':
                kickers.append(el)
        return kickers

