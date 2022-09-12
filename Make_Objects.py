class Objects():
    def __init__(self,file=r"Twiss_2_Alpha_Adjusted_Map_Tune.g4bl",ebtoggle=1,nevents=300,tune=True,zspacing=30,profile="profile.txt",fringe=True,**kwargs): #Setting up the class with things we ALWAYS need to be in the file - Note profile is called by default "profile.txt", you can either change this in the object call or here and use a default, similar for nevents etc - the choice is yours.
        self.fringe = fringe
        self.zspacing = zspacing
        self.file=file
        self.tune_bool=tune
        self.n_events=nevents
        self.profile = profile
        self.ebtoggle=ebtoggle #set to one for equads
        self.load_beam = self.Beam_Setup()
        self.load_generics = self.Generic_Setup()

    def Open_File(self):    #This writes the basic intro to the file
        f = open(self.file,"w+")
        f.write(self.load_beam)
        f.write("\n \n ##BEAM LOADED - LOADING GENERICS## \n")
        f.write(self.load_generics)
        f.close()

    def Read_Objects(self,twiss_objects):   #Reads in and #prints objects from passed from the main method
        for i in twiss_objects:
            print(i['name'],i['k1L'],i['S'],i['type'])

    def Beam_Setup(self):   #Handles the beam setup, can be changed for your desired setup here
        str = "##BEAM SETUP CAN BE HANDLED IN MAKE_OBJECTS ##\nphysics QGSP_BERT disable=Decay \n"
        str1= 'g4ui when=4 "/vis/viewer/set/lineSegmentsPerCircle 80"\ng4ui when=4 "/vis/scene/add/axes -250 0 7200 500 mm"\ng4ui when=4 "/vis/viewer/set/viewpointVector -9.77e-08 1 4e-10"\n'
        str2= 'param M=938.272 KE=0.1 # beam in MeV\nparam P=sqrt(($M+$KE)*($M+$KE)-$M*$M)  # Calculate the momentum.\n'
        str3='param minStep=5\nparam maxStep=10\nparam histoFile=g4beamlineFullAlpha.root\nparam shift1=0 #-2000 #shift all elements by this amount\nparam trackthese=100\n#expandworld 8000,5000,9000\nparam killj=1\n'
        Beam1= "beam ascii filename=/home/steve/Automatch_No_Spread/matched_beam.G4Beam nEvents=%i \n" %(self.n_events) 	#NOTE, this line uses the beam constructed from "Beam_Prep.py" to create events used within G4Beamline - you can instead set up a generic beam however I reccomend doing it this way so the beam stays identical run to run.
        Beam2='reference referenceMomentum=13.699066 particle=anti_proton beamZ=0.0\n'
        final='param shift=0\n'
        strings = str + str1 + str2+ str3 + Beam1 + Beam2 + final
        return strings

    def Generic_Setup(self):    #Handles setting up generic objects, again this needs modifiying, for example focusquad here is modular, with an iron length of 100, field length 100 and radius 60. Naturally your quads will be different, this is where we construct each modular object type! Be careful in this construction, it defines your entire program!
        if self.fringe == True:
            focus_quad = "genericquad focusquad fieldLength=100 ironLength=100 ironRadius=60 \\\napertureRadius=30 ironColor=0,0,1 kill=1 fringe='-0.7233,10.39,1.003,-6.39,0.8186,2.049'\n"
            defocus_quad = "genericquad defocusquad fieldLength=100 ironLength=100 ironRadius=60 \\\napertureRadius=30 ironColor=1,0,0 kill=1 fringe='-0.7233,10.39,1.003,-6.39,0.8186,2.049'\n"
        else:
            focus_quad = "genericquad focusquad fieldLength=100 ironLength=100 ironRadius=60 \\\napertureRadius=30 ironColor=0,0,1 kill=1\n"
            defocus_quad = "genericquad defocusquad fieldLength=100 ironLength=100 ironRadius=60 \\\napertureRadius=30 ironColor=1,0,0 kill=1\n"
        drift = "genericquad drift fieldLength=100 ironLength=100 ironRadius=60 \\\napertureRadius=30 ironColor=0.5,0.5,0 kill=1\n"
        corrector = "genericbend corrector fieldLength=37 fieldWidth=60 fieldHeight=60 \\\nironHeight=70 ironLength=37 ironWidth=70 ironColor=1,0,1 kill=1\n"
        turnedoffquad = "genericquad offquad fieldLength=100 ironLength=100 ironRadius=60 \\\napertureRadius=30 ironColor=0,0.5,0.5 kill=1\n"
        detector= "virtualdetector detector material=Vacuum radius=200.0 color=0,1,0 referenceParticle=1 coordinates=reference format=ascii\n"
        strings = focus_quad + defocus_quad + drift + corrector + turnedoffquad + detector
        return(strings)

    def Place_Objects(self,twiss_objects): #This places all objects sequentially, you can see it calls methods specific to a quad, bend etc. You can add more based on how many "types" you have, for example if you have 2 quadrupole types in the line i suggest making a second function to handle "QUADRUPOLE2".
        f=open(self.file,"a+")
        f.write("\n \n ## GENERICS LOADED - LOADING IN BEAMLINE ELEMENTS ## \n \n")
        t = []
        for i in twiss_objects:
            if i['type'] == 'QUADRUPOLE':
                self.Make_Quadrupole(i,f)
            elif i['type'] == 'HKICKER':
                self.Make_Hkicker(i,f)
            elif i['type'] == 'VKICKER':
                self.Make_Vkicker(i,f)
            elif i['type'] == 'MATRIX':
                self.Make_Bend(i,f)
            elif i['type'] == 'MONITOR':
               self.Make_Monitor(i,f)
            elif i['type'] == "EndLine04":
                self.Make_End(i,f)
            else:
                t.append([i['name'],i['type']])

        f.close()

    def Make_Quadrupole(self,quad,file):	#These are set up currently for a custom G4Beamline build that has gradients defined a Egradient and Bgradient, you will need to change this in your version! It is vitally important you do so!
        if float(quad['grad'])>0:
            if self.ebtoggle == 1:
                str = "place focusquad z=%s rename=%s Egradient=%s rotation=Z45\n" % (quad['S'],quad['name'],quad['grad'])
            else:
                str = "place focusquad z=%s rename=%s Bgradient=%s \n" % (quad['S'],quad['name'],quad['grad'])
        elif float(quad['grad'])<0:
            if self.ebtoggle == 1:
                str = "place defocusquad z=%s rename=%s Egradient=%s rotation=Z45\n" % (quad['S'],quad['name'],quad['grad'])
            else:
                str = "place defocusquad z=%s rename=%s Bgradient=%s \n" % (quad['S'],quad['name'],quad['grad'])
        else:
            if self.ebtoggle == 1:
                str ="place offquad z=%s rename=%s Egradient=%s rotation=Z45\n" % (quad['S'],quad['name'],quad['grad'])
            else:
                str = "place offquad z=%s rename=%s Bgradient=%s \n" % (quad['S'],quad['name'],quad['grad'])
        file.write(str)

    def Make_Hkicker(self,object,file):   #Making Horizontal Kickers

  

        name = object['name']
        name = name.replace(".","_")
        dist = object['S']
        tune_string= "tune %s_cor z0=0 z1=(%s+180) initial=0 step=0.0000001 expr=Py1/Pz1 tolerance=0.001 maxIter=1000\n" %(name,object['S'])
        #tune_string= "tune %s_cor z0=0 z1=(%s+180) initial=0 step=0.0000001 expr=Px1/Pz1 tolerance=0.001 maxIter=1000\n" %(name,object['S'])

        str = "place corrector z=%s rename=%s By=%s_cor*0.000001\n\n"%(object['S'],name,name) #By=%s_cor*0.000001\n
        if self.tune_bool==True:
            string=tune_string+str
        else:
            str = "place corrector z=%s rename=%s \n"%(object['S'],name)
            string=str
        file.write(string)

    def Make_Vkicker(self,object,file):	#Making Vertical kickers

        name = object['name']
        name = name.replace(".","_")
        dist = object['S']
        tune_string= "tune %s_cor z0=0 z1=(%s+180) initial=0 step=0.0000001 expr=Px1/Pz1 tolerance=0.001 maxIter=1000\n" %(name,object['S'])
        if name=="LNE_ZCV_0427":
            tune_string = "tune %s_cor z0=0 z1=(%s+100) initial=0 step=0.0000001 expr=x1 tolerance=0.01 maxIter=1000\n" %(name,object['S'])
        #tune_string= "tune %s_cor z0=0 z1=(%s+180) initial=0 step=0.0000001 expr=Py1/Pz1 tolerance=0.001 maxIter=1000\n" %(name,object['S'])

        str = "place corrector z=%s rename=%s By=%s_cor*0.000001\n\n"%(object['S'],name,name) #By=%s_cor*0.000001\n
        if self.tune_bool==True:
            string=tune_string+str
        else:
            str = "place corrector z=%s rename=%s \n"%(object['S'],name)
            string=str
        file.write(string)


    def Make_Bend(self,object,file):	#The bends used here use externally generated field maps, you can use a generic bend should you wish - this is the part which will require the most changing depending on each user!

        if object['name'] == "LNR.ZDFA.0610":	#Here I specifically call each bending object by its name and construct them as special objects - you will need your own workaround! You could use a keyword call for i.e "MATRIX" type and construct generic bends, thus doing everything in 1 "if" loop - it should be an easy change. Any issues just contact me!
            str1= "box curvebox height=130 width=5 length=400 kill=0\nparam Emag=0.11\nparam angE=3.9465\nparam " \
                  "Exf=$Emag*cos($angE*deg)\nparam Ezf=$Emag*sin($angE*deg) \n"
            str2= "fieldexpr zdfafield Ex=-$Exf Ez=$Ezf length=399.0515 width=150 height=115 nX=100 nY=10 nZ=100\n"
            str3= "fieldmap map610 file=/home/steve/PDF_2_G4BL/New_Lines/zdfa0610_3mm.B\n"
            str4= "place map610 x=0 y=0 z=%s gradient=0.9010/1000000" \
                  "z=0->10600 linux \nlabel text=610 coordinates=%s size=56\n" %(object['S'],object['S'])
            str5= "cornerarc z=%s angle=12.605 centerRadius=1818.19213 #[ (360/12.605)*400mm =circumf, /(2pi) =radius]\n" % (object['S'])
            string = str1 + str2 + str3 + str4 +str5
            file.write(string)

        if object['name'] == "LNE.ZDFHR.0127":
            str1= "box curvebox height=130 width=5 length=400 kill=0\nparam Emag=-0.11\nparam angE=3.9465\nparam " \
                  "Exf=$Emag*cos($angE*deg)\nparam Ezf=$Emag*sin($angE*deg) \n"
            str2= "fieldexpr zdfafield Ex=-$Exf Ez=$Ezf length=399.0515 width=150 height=115 nX=100 nY=10 nZ=100\n"
            str3= "fieldmap map127 file=/home/steve/PDF_2_G4BL/New_Lines/zdfhr0127_3mm_2.B\n"
            str4= "place map127 x=0 y=0 z=%s gradient=0.8746/1000000" \
                  "z=0->10600 linux \n" %(object['S'])
            str5= "cornerarc z=%s angle=-12.605 centerRadius=1818.19213 #[ (360/12.605)*400mm =circumf, /(2pi) =radius]\n" %(object['S'])
            string = str1 + str2 + str3 + str4 +str5
            file.write(string)

        if object['name'] == "LNE.ZDSHR.0132":
            str1= "param len0132=(33.16*(pi/180)*600)\nparam lenh0132=173.625\n"    #calculating the arc length
            str2= "fieldmap map132 file=/home/steve/PDF_2_G4BL/New_Lines/zdshr0132_3mm_2.B\n"    #Loading map with adjusted graident
            str3= "place map132 z=%s gradient=0.9173/1000000\n" %(object['S'])  #Placing said map
            str4= "cornerarc z=%s centerRadius=600 angle=-33.16\n" %(object['S'])    #This may need to go negative not positive addition? Corner arc also adjusts centreline coords
            string=str1+str2+str3+str4
            file.write(string)

        if object['name'] == "LNE.ZDSHR.0415":
            name = object['name']
            name = name.replace(".", "_")
            str1= "param len0415=50.42*(pi/180)*600\nparam lenh0415=($len0415/2)\n"
            str2= "fieldmap map415 file=/home/steve/PDF_2_G4BL/New_Lines/zdshr0415_3mm_2.B\n"
            tunestring = "tune %s_MAPTUNE z0=0 z1=33114.5 initial=0.9 step=0.00001 expr=x1 tolerance=0.000001 maxIter=10000\n" % (name)
            str3= "place map415 z=%s-251.85 gradient=0.000001*%s_MAPTUNE\n" %(object['S'],name)

            str4= "cornerarc z=%s-251.85 centerRadius=600 angle=-50.42\n" %(object['S'])
            string=str1+str2+tunestring+str3+str4
            file.write(string)

        if object['name'] == "LNE.ZDSHR.0220":
            tempdet = "place detector z=%s-251.85 rename=tempdet \n" %(object['S'])
            name = object['name']
            name = name.replace(".", "_")
            str1= "param len0220=45.77*(pi/180)*600\nparam lenh0220=($len0220/2)\n"
            str2= "fieldmap map220 file=/home/steve/PDF_2_G4BL/New_Lines/zdshr0415_3mm_2.B\n"    #Using field map from 415 as it is the same component
            tunestring = "tune %s_MAPTUNE z0=0 z1=24552.45483 initial=0.9 step=0.00001 expr=x1 tolerance=0.000001 maxIter=10000\n" % (name)
            str3= "place map220 z=%s-251.85 rotation=Y0.8 gradient=0.000001*%s_MAPTUNE\n" %(object['S'],name)
            str4= "cornerarc z=%s-251.85 centerRadius=600 angle=-45.77\n" %(object['S'])
            string=tempdet+str1+str2+tunestring+str3+str4
            file.write(string)


    def Make_Drift(self,object,file):	#Inserts drifts
        file.write("#Temp Drift \n")
    def Make_Monitor(self,object,file): #Places virtual derectors
        string="place detector z=%s rename=%s \n" %(object['S'],object['name'])
        file.write(string)
    def Make_End(self,object,file):	#Closes the beam - produces a profile used for analysis.
        string="place detector z=%s rename=%s \n" %(object['S'],"Handover") #object['type']
        profile = "profile file=%s particle=anti_proton zloop=0:%s:%s coordinates=centreline\n" % (self.profile,object['S'],self.zspacing)
        file.write(string + profile)


