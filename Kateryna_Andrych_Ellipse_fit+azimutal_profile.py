import matplotlib.pyplot as plt
import numpy as np
import fnmatch
import os
from astropy.io import fits
from skimage.measure import EllipseModel
from matplotlib.patches import Ellipse
from scipy import interpolate
import matplotlib.patches as patches
import functions as f

"""
Downloads fit file as numpy array, fits ellipse to pixels with brightness >=brightnesslimit and plot brightness profile along the fitting using bilinear interpolation for each point.

"""

dirdat = '/media/kateryna/Data_Lin/PhD/IRDAP_reduction/kluska_0102.D-0696(A)/0-3_reduction/'  #where data
fig_folder= '/media/kateryna/Data_Lin/PhD/IRDAP_reduction/kluska_0102.D-0696(A)/0-3_reduction/' #where save images
stars = ['iras08544-4431_calib','hr4049','iras15469-5311','iras17038-4815','iw_car_calib','ru_cen_calib','u_mon_2019-01-03_calib','u_mon_2019-01-14_calib','u_mon_combined'] #objects names that are used as data folder name 

fittypes=['Q_phi'] 
annuluses=['0-3']
brightnesslimit=15

    

for star in stars:
    print(star)
    for annulus in annuluses:
        print('Annulus for stellar polarisation '+annulus)
        savefig=fig_folder+star+'/analysis/'
        try:
        # Create target Directory
           os.mkdir(savefig)
        except FileExistsError:
           print("Directory " , savefig ,  " already exists")

        savefig=savefig+'azimutal_profile/'
        try:
        # Create target Directory
           os.mkdir(savefig)
        except FileExistsError:
           print("Directory " , savefig ,  " already exists")

        
        logfile=open(savefig+star+'_'+annulus+'_ellipse_logfile.txt','w+')

        dirName=savefig+'pixels_for_interpolation/'
        try:
        # Create target Directory
           os.mkdir(dirName)
        except FileExistsError:
           print("Directory " , dirName ,  " already exists")
    
        savefig2=dirdat+'all_analysis_results/'
        try:
        # Create target Directory
           os.mkdir(savefig2)
        except FileExistsError:
           print("Directory " , savefig2 ,  " already exists")

        savefig2=savefig2+'azimutal_profile/'
        try:
        # Create target Directory
           os.mkdir(savefig2)
        except FileExistsError:
           print("Directory " , savefig2 ,  " already exists")

        
        for fittype in fittypes:
            image,hdul= f.Loadimages(star,fittype,annulus,dirdat)
            n = image.shape[0]    
            ps = 12.27 #mas per pixel for IRDIS
            
            
            #Creating grid for the multiplying the image by the separation from star.
            xr = np.linspace(-n/2, n/2, num=n)
            yr = np.linspace(-n/2, n/2, num=n)
            x0 = 0.5
            y0 = 0.5
            xr = xr-x0
            yr = yr-y0
            Xr, Yr = np.meshgrid(xr, yr)
            R = np.sqrt(Xr**2 + Yr**2)
            
         
            image_an= image*1.0
            n = image_an.shape[0]
            

            image_an = image_an*(R<30)
        
            #selecting points that will be used for the ellipse fitting 
            data_p= np.where(image_an >= brightnesslimit)
            data_p1=np.vstack((data_p[1],data_p[0])).transpose()    
            a_points = data_p1
            x = a_points[:, 0]
            y = a_points[:, 1]

            #fiting an ellipse
            ell = EllipseModel()
            ell.estimate(a_points)
            xc, yc, a, b, theta = ell.params
            shift=n/2-0.5
            #calculating points of ellipse based on the previously defined parameters
            points=ell.predict_xy(np.linspace(0, 2 * np.pi, 50),params=(xc,yc,a,b,theta)) #without scaling into mas, only pixels` numbers
            points_mas=ell.predict_xy(np.linspace(0, 2 * np.pi, 50),params=((xc-shift)*ps,(yc-shift)*ps,a*ps,b*ps,theta)) #in mas
            
            sigma_2=np.sum(ell.residuals(a_points)**2)
            #writting results into the logfile
            logfile.writelines(star+'\n')
            logfile.writelines('Annulus for stellar polarisation '+annulus+'\n')
            logfile.writelines("center = (%f , %f) \n" % (xc, yc))
            logfile.writelines("angle of rotation = %f \n" % theta)
            logfile.writelines("half axes im mas= %f, %f \n" % (a*ps,b*ps))
            logfile.writelines("sigma^2 = %f \n" % sigma_2)
            
            #plotting only points and ellipse
            fig, axs = plt.subplots(1, 1, sharex=True, sharey=True)
            axs.scatter(x, y)
            axs.plot(xc, yc, "+", color='red')
            plt.plot(points[:,0],points[:,1], color='g')
            plt.title(star+' '+fittype)
            plt.savefig(savefig+star+'_'+annulus+'_el.jpeg',bbox_inches='tight', pad_inches=0.1)
            plt.close()
            
            #image in sinh scale and mas on axes, fitted ellipse, star position pointed by green +, ellipse center - red +  
            lim=10

            image=np.arcsinh(image_an)

            fig, ax = plt.subplots()
            max = np.max(image)
            min=np.min(image)
            d = n * ps / 2
            plt.imshow(image, vmin=min, vmax=max, extent=(-d, d, d, -d))
            plt.plot(0, 0, "+", color="white")
            plt.plot(points_mas[:,0], points_mas[:,1],color='red')
            ax.plot((xc-shift)*ps, (yc-shift)*ps, "+", color='red')
            plt.xlim(-lim * ps, lim * ps)
            plt.ylim(-lim * ps, lim * ps)
            plt.xlabel('mas')
            plt.ylabel("mas")
            plt.colorbar()
            plt.tight_layout      
            plt.title(star+' '+fittype)
            plt.savefig(savefig+star+'_'+annulus+'_im+el.jpeg',bbox_inches='tight', pad_inches=0.1)   
            plt.title(star) 
            plt.close()
            

            #testing
            if a>b:
                pointmajax=ell.predict_xy(np.linspace(0, np.pi, 2),params=(xc,yc,a,b,theta)) #without scaling into mas, only pixels` numbers
                startpositionell=np.linspace(np.pi/2, np.pi*3/2, 2)[pointmajax[:,0]==np.max(pointmajax[:,0])]
                arrowendarray=ell.predict_xy(np.linspace(0+np.deg2rad(30), np.pi+np.deg2rad(30), 2),params=(xc,yc,a,b,theta))
                arrowend=arrowendarray[pointmajax[:,0]==np.max(pointmajax[:,0]),:]
                print(arrowend)
                pointstart=pointmajax[pointmajax[:,0]==np.max(pointmajax[:,0]),:]
                print(pointstart)
            else:
                pointmajax=ell.predict_xy(np.linspace(np.pi/2, np.pi*3/2, 2),params=(xc,yc,a,b,theta)) #without scaling into mas, only pixels` numbers
                startpositionell=np.linspace(np.pi/2, np.pi*3/2, 2)[pointmajax[:,0]==np.max(pointmajax[:,0])]
                arrowendarray=ell.predict_xy(np.linspace(np.pi/2+np.deg2rad(30), np.pi*3/2+np.deg2rad(30), 2),params=(xc,yc,a,b,theta))
                arrowend=arrowendarray[pointmajax[:,0]==np.max(pointmajax[:,0]),:]
                print(arrowend)
                pointstart=pointmajax[pointmajax[:,0]==np.max(pointmajax[:,0]),:]
                print(pointstart)
            
            startangle=np.arctan2((pointstart[0,1]-yc), (pointstart[0,0]-xc))
            print(startangle)
            #plotting brightness profile along the fitting line   

            #position angle is mesured from the east counterclockwise   
            
            #selecting pixels that ellipse goes through         
            mask = np.zeros((n,n))
            x = np.int32(np.round(points[:,0]))
            y = np.int32(np.round(points[:,1]))
            mask[y,x] = 255
            ellprofile=image_an[mask==255]
            data_el= np.where(mask == 255)

            #coordinates of this pixels and position angle for each of them
            el_p=np.vstack((data_el[1],data_el[0])).transpose()   
            x_el = el_p[:, 0]
            y_el = el_p[:, 1]
            r, pos_angle = f.cart2polar(x_el, y_el,xc,yc,startangle)
            pos_angle=np.rad2deg(pos_angle)
            #print('look')
            #print(pos_angle)
            indices=pos_angle.argsort()
            #print('sorted')
            #print(pos_angle[indices])
            #print((el_p[indices,:]-shift)*ps)
            #image in sinh scale and mas on axes, fitted ellipse, black dots point pixels that ellipse goes through 
            lim=10




            #plotting brighness profile along the ellipse using original pixel values
            plt.plot(pos_angle[indices],ellprofile[indices], marker='d')
            plt.title(star)
            plt.xlabel('position angle, deg')
            plt.ylabel("arbitrary counts")
            plt.savefig(savefig+star+'_'+annulus+'_profile_per_fit.jpeg',bbox_inches='tight', pad_inches=0.1)
            plt.close()




            #creating brifhtness profile along the ellipse using bilinear interpolation and regulary spaced points
            npoints=30
            #if a>b:
            ellangle=np.linspace(0, 2 * np.pi, npoints)-startangle
            
            for i in range(0,len(ellangle)):
                if ellangle[i]<0:
                    ellangle[i]=2*np.pi+ellangle[i]
            points_ellipse=ell.predict_xy(np.linspace(0, 2 * np.pi, npoints),params=(xc,yc,a,b,theta))
            points_ellipse_plot=ell.predict_xy(np.linspace(0, 2 * np.pi, npoints),params=((xc-shift)*ps,(yc-shift)*ps,a*ps,b*ps,theta)) #in mas
            inter_prof=np.zeros(points_ellipse.shape[0])
            ii=0
            inter_pos_angle=np.zeros(points_ellipse.shape[0])
            for xcoor,ycoor in points_ellipse:
                xint=int(xcoor-xcoor%1) 
                yint=int(ycoor-ycoor%1)
                xshift=1
                yshift=1
                
                n_pix=[(xint,yint,image_an[xint,yint]),(xint,yint+yshift,image_an[xint,yint+yshift]),(xint+xshift,yint,image_an[xint+xshift,yint]),(xint+xshift,yint+yshift,image_an[xint+xshift,yint+yshift])]
              
                ncoordx=np.array([xint,xint,xint+xshift,xint+xshift])
                ncoordy=np.array([yint,yint+yshift,yint,yint+yshift])
                gifimage=plt.imshow(image, vmin=min, vmax=max,extent=(-d, d, d, -d))
                plt.plot(points_ellipse_plot[:,0], points_ellipse_plot[:,1],color='red')
                plt.xlim(-lim * ps, lim * ps)
                plt.ylim(-lim * ps, lim * ps)

                plt.colorbar()
                plt.scatter((ncoordx-shift)*ps,(ncoordy-shift)*ps,color='black',s=20)
                plt.scatter((xcoor-shift)*ps,(ycoor-shift)*ps,color='red', s=40)
                pospoint = np.arctan2((ycoor-yc), (xcoor-xc))-startangle
                if pospoint<0:
                    pospoint=2*np.pi+pospoint
                plt.title(star+str(np.rad2deg(pospoint)))
                plt.xlabel('mas')
                plt.ylabel("mas")
                if ii<10: plt.savefig(dirName+star+'_'+annulus+'_profile_pixels0'+str(ii)+'.jpeg',bbox_inches='tight', pad_inches=0.1)
                else:plt.savefig(dirName+star+'_'+annulus+'_profile_pixels'+str(ii)+'.jpeg',bbox_inches='tight', pad_inches=0.1)


                plt.close()
               
                
                inter_prof[ii]=f.bilinear_interpolation(xcoor,ycoor,n_pix)
                inter_pos_angle[ii]=pospoint
                ii+=1
            

            
            inter_pos_angle=np.rad2deg(inter_pos_angle)
            inter_indices=inter_pos_angle.argsort()
            print(inter_pos_angle)
                  
            plt.plot(inter_pos_angle[inter_indices],inter_prof[inter_indices], 'ro-',label='interpolated')
            plt.plot(pos_angle[indices],ellprofile[indices],'bo:',label='pixel_values')
            plt.xlabel('position angle, deg')
            plt.ylabel("arbitrary counts")
            plt.title('Profile along the ellipse')
            plt.legend()
            plt.savefig(savefig+star+'_'+annulus+'_profile_compare.jpeg',bbox_inches='tight', pad_inches=0.1)
            plt.savefig(savefig2+star+'_'+annulus+'_profile_compare.jpeg',bbox_inches='tight', pad_inches=0.1)
            plt.close()

            plt.plot(inter_pos_angle[inter_indices],inter_prof[inter_indices], 'ro-',label='interpolated')
            plt.xlabel('position angle, deg')
            plt.ylabel("arbitrary counts")
            plt.title('Profile along the ellipse')
            plt.legend()
            plt.savefig(savefig+star+'_'+annulus+'_bilinearprofile.jpeg',bbox_inches='tight', pad_inches=0.1)
            plt.close()
            logfile.writelines('\n')
            logfile.close()
            
            d = n * ps / 2
            plt.imshow(image, vmin=min, vmax=max, extent=(-d, d, d, -d))
            plt.plot(points_mas[:,0], points_mas[:,1],color='red')
            plt.xlim(-lim * ps, lim * ps)
            plt.ylim(-lim * ps, lim * ps)
            plt.colorbar()
            plt.plot((pointstart[:,0]-shift)*ps,(pointstart[:,1]-shift)*ps,color='red',marker="o",markersize=4)
            #plt.scatter((x_el[indices[0]]-shift)*ps,(y_el[indices[0]]-shift)*ps,color='blue')
            style = "Simple, tail_width=0.5, head_width=8, head_length=8"
            kw = dict(arrowstyle=style, color="red")
            arrowxend,arrowyend=arrowend[:,0],arrowend[:,1]
            #print(arrowxend[0],arrowyend[0])
            arrowxstart,arrowystart=pointstart[:,0],pointstart[:,1]
            #print(arrowxstart[0],arrowystart[0])
            arrow=patches.FancyArrowPatch(((arrowxstart[0]-shift)*ps,(arrowystart[0]-shift)*ps), ((arrowxend[0]-shift)*ps,(arrowyend[0]-shift)*ps),connectionstyle="arc3,rad=0.2", **kw)     
            plt.gca().add_patch(arrow)
            #plt.scatter((x_el-shift)*ps,(y_el-shift)*ps,color='green')      
            plt.title(star)
            plt.xlabel('mas')
            plt.ylabel("mas")        
            plt.savefig(savefig+star+'_'+annulus+'_profile_pixels.jpeg',bbox_inches='tight', pad_inches=0.1)
            plt.savefig(savefig2+star+'_'+annulus+'_profile_pixels.jpeg',bbox_inches='tight', pad_inches=0.1)
            plt.close()
      




            
           
            

           
