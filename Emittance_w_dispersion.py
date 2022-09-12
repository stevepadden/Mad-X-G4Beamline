#As the automatic construction routines produce many BPM files, when the G4Beamline code runs it populates those files with data. This script runs through them all and plots them, it is very useful for visualation of how the beam travels.

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sea

# Note to the user, this program -MAY- throw an error when saving the files for each BPM on windows, there exist ...
# solutions online regarding changing folder types, however i reccomend not bothering with this, simply delete any ...
# files begining with .png in the folder specific for the BPMS, this then should have no issues and removes ...
# the need to worry about setting folder permissions to numeric (not simple). On linux this should not be an ...
# issue - I point blank refuse to use a Mac to test due to chubby fingers not liking the command button. Regards \S

#dir = "/home/steve/BPMS"  # File directory containing the BPM files
dir = "/home/steve/Automatch_No_Spread" #Adjust to your directory where the G4BL bpm files are!
dats = ["x", "y"]# Flick between x or y at convinience - does need formatting as a text string.
shift_data = False # Centers the ellipse back on 0,0 for the purposes of plotting, depending on how G4BL is used ...
#Sometimes the data can shift away from 0,0 despite it being centered on 0,0 from the beamline perspective - be careful.
show_plots = True  # Boolean for displaying plots or not.
calc_emit = True
mm_conv = True #Flag to convert to mm from m, set true to use mm
avg_momentum = 13.699065661570090 #Adjust to your beam's momentum here, necessary to ensure correct calculation of emittance.
df_disp = pd.read_pickle("/home/steve/Automatch_No_Spread/dispersion.pkl")	#Reads the dispersion calculated via mad-x


ermsx = []
ermsy =[]
sx = []
sy = []
#Calulating emittance in x and y 
def emit(data,flag):
    if flag == "x":

        sxd = data["x"]
        sxprimed = data["xprime"]
        sxxprimed = data["xxprime"]
        x2 = np.power(sxd, 2)
        xprime2 = np.power(sxprimed, 2)
        avgx2 = np.average(x2)
        avgxprime2 = np.average(xprime2)

        avgxxprime = np.average(sxxprimed)
        avgxxprime2 = np.power(avgxxprime, 2)
        emitx = np.sqrt((avgx2 * avgxprime2) - avgxxprime2)
        data["emitx"] = emitx
        data["fname"] = filename
        return (emitx)
    elif flag== "y":

        syd = data["y"]
        syprimed = data["yprime"]
        syyprimed = data["yyprime"]
        y2 = np.power(syd, 2)
        yprime2 = np.power(syprimed, 2)
        avgy2 = np.average(y2)
        avgyprime2 = np.average(yprime2)
        # here we average then square
        avgyyprime = np.average(syyprimed)
        avgyyprime2 = np.power(avgyyprime, 2)
        emity = np.sqrt((avgy2 * avgyprime2) - avgyyprime2)
        data["emity"] = emity
        data["fname"] = filename
        return emity

def findmax(x,y):
    x1 = max(x)
    x2 = abs(min(x))
    y1 = max(y)
    y2 = abs(min(y))
    if x2>x1:
        maxx = x2
        minx = -x2
    else:
        maxx = x1
        minx = -x1
        
    if y2>y1:
        maxy = y2
        miny = -y2
    else:
        maxy = y1
        miny = -y1
    endlim=(max(maxx,maxy))
    return minx,maxx,miny,maxy,endlim
        
    

for filename in os.listdir(dir):
    if "LNE" not in filename:
        if "Handover" not in filename:
            continue
        #if "EndLine" not in filename:
        #    continue
        #if "Handover" not in filename:
        #    continue
    if filename.endswith(".txt"):
        File = os.path.join(dir, filename)  # Setting the file to this iteration
        filename = filename.strip(
            ".txt")  # Removing the .txt for later formatting - if your file isnt a .txt change to the extension here.
        pd.set_option('display.max_columns', None)  # set display columns pandas
        cols = ["x", "y", "z", "Px", "Py", "Pz", "t", "PDGid", "EventID", "TrackID", "ParentID", "Weight"]
        data = pd.read_csv(File, sep=" ", skiprows=3, header=None, names=cols)  # read data with above vols

        #PLOTTING X Y HEATMAP
        cmap ="viridis"
        fig,ax = plt.subplots()
        heatmap, xedges, yedges = np.histogram2d(data["x"], data["y"], bins=50)
        rmsx = np.sqrt(np.sum(np.power(data["x"],2) )/len(data["x"]) )
        print("RMSX")
        print(rmsx)
        rmsy = np.sqrt(np.sum( np.power (data["y"] ,2) ) / len(data["y"]))
        print("RMSY")
        print(rmsy)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        ax.set_title(filename)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        minx,maxx,miny,maxy,endlim = findmax(data["x"],data["y"])

        ax.set_xlim((-endlim,endlim))
        ax.set_ylim((-endlim,endlim))
        tt=ax.imshow(heatmap.T,extent=extent,origin='lower',cmap=cmap,aspect="equal")
        ax.set_facecolor('#440154')
        str = "$\sigma_{RMSX} =$ %.2f (mm) \n$\sigma_{RMSY} =$ %.2f (mm)" %(rmsx,rmsy)
        ax.plot([], [], ' ', label=str)
        ax.legend(frameon=False,facecolor = '#440154', labelcolor="w")
        cbar=plt.colorbar(tt,ax=ax)

        if show_plots == True:
            fig.show()
        fig.savefig(dir+"/"+filename+"Spatial_2hist.png")
        # PLOTTING PX PY HEATMAP
        fig,ax = plt.subplots()
        heatmap, xedges, yedges = np.histogram2d(data["Px"], data["Py"], bins=50)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        ax.set_title(filename)
        ax.set_xlabel("Px (MeV/c)")
        ax.set_ylabel("Py (MeV/c)")
        minx,maxx,miny,maxy,endlim = findmax(data["Px"],data["Py"])

        ax.set_xlim((-endlim,endlim))
        ax.set_ylim((-endlim,endlim))
        tt=ax.imshow(heatmap.T, extent=extent, origin='lower',cmap=cmap,aspect="equal")
        ax.set_facecolor('#440154')
        fig.colorbar(tt,ax=ax)

        if show_plots == True:
            fig.show()
        fig.savefig(dir+"/"+filename+"Momentum_2hist.png")

        data["x"]=data["x"]/1000
        data["y"]=data["y"]/1000
        data["z"]=data["z"]/1000

        z = data["z"]
        z0 = z[0]

        idx = df_disp['s'].sub(z0).abs().idxmin()
        temp = df_disp.loc[[idx]]

        totp = np.sqrt(data["Pz"]**2 + data["Px"]**2 + data["Py"]**2)
        delp = (data["Pz"]- avg_momentum) / avg_momentum

        dispersion = float(temp["DX"])
        dispersionprime = float(temp["DDX"])


        data["x"] = data["x"]#+(delp*dispersion)
        xprime = ((data["Px"]) / data["Pz"])#+(delp*dispersionprime)# Calculating x' using for small angles (x1-x2/L) ~ (Px/Pz) assuming paraxial approximation
     

        data["xprime"] = xprime
        yprime = data["Py"] / data["Pz"]
        data["yprime"] = yprime

        if mm_conv == True:
            data["xprime"] = data["xprime"] *1000
            data["yprime"] = data["yprime"] *1000
            data["x"] = data["x"] *1000
            data["y"] = data["y"] *1000

        data["xxprime"] = data["x"] * data["xprime"]  #- (delp*dispersionprime)
        data["yyprime"] = data["y"] * data["yprime"]
        if shift_data == True: #Shifting the data to 0,0 if needed

            avgx = data["x"].mean()  # X can sometimes be shifted, global coords - just moving it back to roughly 0 for the purposes of plotting. Draw attention to this shift if required
            data["x"] = data["x"] - avgx
            avgy = data["y"].mean()
            data["y"] = data["y"] - avgy


        import matplotlib.pyplot as plt
        for xory in dats:
            from matplotlib.patches import Ellipse
            x = data[xory]  # setting our "x" data to "y"
            #y= data["P"+xory]
            y = data[xory + "prime"]  # and y data to y'
            xmax = max(abs(x)) * 1.5#1.012
            ymax = max(abs(y)) * 1.5#1.012
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            plt.xlim(-xmax, xmax)
            plt.ylim(-ymax, ymax)

            ax.spines['left'].set_position(
                'zero')  # ax setting to center the axis, comment any ax.spines to return to edge axis.
            ax.spines['bottom'].set_position('zero')

            # Eliminate upper and right axes
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')

            # Show ticks in the left and lower axes only
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.set_xlabel(xory, fontsize=15)  # Setting the labels for x and y
            ax.set_ylabel(xory + "'", rotation=0, fontsize=15)
            ax.xaxis.set_label_coords(0.98, 0.45)  # the X label position
            ax.yaxis.set_label_coords(0.55,
                                      0.95)  # And the y label position (0,0 bottom left, 0.5,0.5 center, 1,1 top right)
            def eigso(cov):  # Function to return eigenvalues of covariant matrix
                vs, vc = np.linalg.eigh(cov)
                ord = vs.argsort()[::-1]
                return vs[ord], vc[:, ord]
            plt.scatter(x, y, marker="+",color="royalblue")  # scatter before ellipse, makes it clearer
            stdevs = [1, 2, 3]  # Standard deviations
            devs = ["68%", "95%", "99.7%"]  # Each stdev amount
            colours = ["red", "black", "orange"]  # colours of ellipses
            for i in range(len(stdevs)):  # calculating and plotting each ellipse.
                ax = plt.subplot(111)
                cov = np.cov(x, y)
                vs, vc = eigso(cov)
                t = np.degrees(np.arctan2(*vc[:, 0][::-1]))
                w, h = 2 * stdevs[i] * np.sqrt(vs)
                if calc_emit == True:
                    e = w * h
                    str22 = devs[i] +(" ($\epsilon_%s = %1.2e$)" % (xory,e))
                    #print(str22)
                else:
                    str22 = devs[i]

                ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                              width=w, height=h,
                              angle=t, color=colours[i], label=str22)
                ell.set_facecolor('none')

                ax.add_patch(ell)
            if calc_emit == True:
                erms = emit(data,xory)
                plt.legend(loc="best",title='$\epsilon_{%s rms} = %.2e$' % (xory, erms),frameon=False)
                if xory == "x":
                    ermsx.append(erms)
                    sx.append(z[0])
                if xory == "y":
                    ermsy.append(erms)
                    sy.append(z[0])

                #plt.legend()
            else:
                ax.legend()
            plt.title(filename)
            plt.savefig(dir +"/" + filename + "_emittance_" + xory + "DISP.png")
            #print(dir +"/" + filename + "_emittance_" + xory + ".png")
            if show_plots == True:
                plt.show()
            plt.clf()



sx,ermsx = zip(*sorted(zip(sx,ermsx)))
sy,ermsy = zip(*sorted(zip(sy,ermsy)))
print(ermsx[0])
print(ermsy[0])
ypos = 3.25

plt.rcParams.update({'font.size':15})
fig,ax = plt.subplots(figsize=(10,7))
ax.plot(sx,ermsx,color="red",label="$\epsilon_{xrms}$")
ax.plot(sy,ermsy,color="blue",label="$\epsilon_{yrms}$")
ax.axvline(0.6,color="black")
ax.text(0.9,ypos,'ZDFA.0610',rotation=90)
ax.axvline(17.498,color="black")
ax.text(16.5,ypos,'ZDFHR.0127',rotation=90)
ax.axvline(18.058,color="black")
ax.text(18.358,ypos,'ZDSHR.0132',rotation=90)
ax.axvline(31.782,color="black")
ax.text(32.082,ypos,'ZDSHR.0415',rotation=90)
ax.axhline(ermsx[0],linestyle="--",color="red",alpha=0.4)
ax.axhline(ermsy[0],linestyle="--",color="blue",alpha=0.4)
ax.legend(frameon=False)
ax.set_title("Statistical emittance along curvilinear beamline axis")
#ax.set_ylim([2.5,5])

ax.set_xlabel("Curvilinear Distance (m)")
ax.set_ylabel("$\epsilon_{rms}$")
fig.show()


fig,ax = plt.subplots(figsize=(10,7))
ax.plot(df_disp["s"],df_disp["DX"],color="green")
ax.axhline(0,color="red")
ax.set_xlabel("Curvilinear Distance , s (m)")
ax.set_ylabel("DX (m)")
ax.set_title("Dispersion from ELENA to ALPHA")
ax.axvline(0.6,color="black")
ax.text(0.9,6,'ZDFA.0610',rotation=90)
ax.axvline(17.498,color="black")
ax.text(16.5,6,'ZDFHR.0127',rotation=90)
ax.axvline(18.058,color="black")
ax.text(18.358,6,'ZDSHR.0132',rotation=90)
ax.axvline(31.782,color="black")
ax.text(32.082,6,'ZDSHR.0415',rotation=90)
plt.show()
