import Tkinter, tkFileDialog, os, ttk, tkMessageBox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib as mpl
import numpy as np
import re
from shirley import shirley
import gaussian
import savitzky
from scipy import interpolate
import colormaps

class MainWindow(Tkinter.Frame):
    def __init__(self,parent=None):
#Creating frames
        Tkinter.Frame.__init__(self)
        self.pack(expand=True,fill='both')
        self.menu_frame = Tkinter.Frame(self)
        self.menu_frame.pack(side=Tkinter.TOP)
        self.image_frame = Tkinter.Frame(self)
	self.image_frame.pack(side=Tkinter.TOP, expand=True, fill='both')
        self.inter_frame = Tkinter.Frame(self)
        self.inter_frame.pack(expand=True)
	self.info_frame = Tkinter.Frame(self)
        self.info_frame.pack(side=Tkinter.TOP, expand=True, fill='both')
        self.slider_frame = Tkinter.Frame(self)
        self.slider_frame.pack(side=Tkinter.TOP)
        self.progress_frame = Tkinter.Frame(self)
        self.progress_frame.pack(side=Tkinter.TOP)
#Change Tkinter themes using ttk
        self.style = ttk.Style()
        self.style.configure('TButton',relief='flat')
        self.available_themes = self.style.theme_names()
#Creating widgets
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig,master=self.image_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side='bottom', fill=Tkinter.BOTH, expand=True)
        self.label = Tkinter.Label(self.info_frame)
        self.label.pack(side=Tkinter.LEFT)
        self.evlabel = Tkinter.Label(self.info_frame)
        self.evlabel.pack(side=Tkinter.RIGHT)
        self.slider = Tkinter.Scale(self.slider_frame,orient=Tkinter.HORIZONTAL,length=300,
                                    command=self.slide_image,to=0,sliderrelief='flat')
        self.slider.pack(side=Tkinter.BOTTOM)
        self.progress = ttk.Progressbar(self.progress_frame,orient="horizontal",length=300, mode="determinate")
        self.progress.pack()

#Creates main application buttons
        self.openbutton = ttk.Button(self.menu_frame, text = "Open", command=self.openfile).pack(side=Tkinter.LEFT)
        self.prevbutton = ttk.Button(self.inter_frame, text = "<", command=self.previmg, state=Tkinter.DISABLED)
        self.prevbutton.pack(side=Tkinter.LEFT)
        self.nextbutton = ttk.Button(self.inter_frame, text = ">", command=self.nextimg, state=Tkinter.DISABLED)
        self.nextbutton.pack(side=Tkinter.LEFT)
        self.testbutton = ttk.Button(self.menu_frame, text = "Double Peak",
                                    command=self.double_peak, state = Tkinter.DISABLED)
        self.testbutton.pack(side=Tkinter.LEFT)
        self.incbutton = ttk.Button(self.menu_frame, text = "+", command=self.increase_psize, state=Tkinter.DISABLED)
        self.incbutton.pack(side=Tkinter.LEFT)
        self.decbutton = ttk.Button(self.menu_frame, text = "-", command=self.decrease_psize, state=Tkinter.DISABLED)
        self.decbutton.pack(side=Tkinter.LEFT)
        self.runbutton = ttk.Button(self.menu_frame, text = "Run", command=self.run, state=Tkinter.DISABLED)
        self.runbutton.pack(side=Tkinter.LEFT)
        self.pizzabutton = ttk.Button(self.menu_frame, text = "Pizza!", command=self.makepizza, state=Tkinter.NORMAL)
        self.pizzabutton.pack(side=Tkinter.LEFT)
        self.configbutton = ttk.Button(self.menu_frame, text = "Configuration",
                                    command=self.config).pack(side=Tkinter.LEFT)
        self.stylebutton = ttk.Button(self.menu_frame, text = "Theme",
                                    command = self.change_style).pack(side=Tkinter.LEFT)
#Canvas event handling and other useful variables
        self.regex = "\d{2,3}(?=\.)"
        self.peak_trigger = False
        self.run_trigger = False
        self.pizza_trigger = False
        self.peak_size = 6
        self.degree = 0
        self.colormap = colormaps.afmcold
        self.toggleViewState = Tkinter.IntVar()
        self.toggleViewState.set(1)
        self.toggleShirleyState = Tkinter.IntVar()
        self.toggleShirleyState.set(1)
        self.toggleGaussianState = Tkinter.IntVar()
        self.toggleGaussianState.set(1)
        self.simmetryState = Tkinter.IntVar()
        self.simmetryState.set(0)
#Define options for opening files
        self.file_opt = options = {}
        options['title'] = "Load your files"
        options['multiple'] = 1
        options['parent'] = self
        self.file_opt2 = options2 = {}
        self.files2 = []
        options2['title'] = "Load your output files"
        options2['parent'] = self
        options2['filetypes'] = [("easyPizza Output Files",("*.xpd","*.hol","*.dat"))]
        options2['initialdir'] = os.getcwd()
#State of the slide bar
        self.Scale = False
#Define the functions called
    def openfile(self):
        """Create a subplot and set index variable to 0 (first file).
        Open multiple text files with the data collected on LNLS, already splited.
        Then, activate the buttons that were Disabled and configure the Slider length."""
        self.ax = self.fig.add_subplot(111)
        self.index = 0
        self.files =  tkFileDialog.askopenfilenames(**self.file_opt)
        if not self.files:
            return
        else:
            self.plotsimple(self.files[0])
            self.Scale = True
            self.evlabel["text"] = os.path.basename(self.files[0])
            self.slider.configure(from_=0,to=len(self.files)-1)
            self.testbutton["state"] = Tkinter.NORMAL
            self.incbutton["state"] = Tkinter.NORMAL
            self.decbutton["state"] = Tkinter.NORMAL
            self.runbutton["state"] = Tkinter.NORMAL
            if len(self.files) > 1:
                self.nextbutton["state"] = Tkinter.NORMAL
                self.prevbutton["state"] = Tkinter.NORMAL

    def plotsimple(self,txt):
        """The axis are cleared, in case any plot was there before.
        The regex searches for the energy in the text file, although this may be unnecessary.
        First the data is smoothed using a Savitzky Golay Filter, to prevent the Shirley to go nuts.
        There are options to use the Shirley Background Subtraction method, and adjust a Gaussian Peak.
        If data has two peaks, there is a Double Peak and + - buttons to adjust the size of the two peaks
        If the Run button is pressed, it integrates the Gaussian Peak adjusted and saves the data on output file."""
        self.ax.clear()
        self.ax.get_xaxis().set_ticks([])
        File = open(txt)
        L = re.findall(r"[-+]?\d*\.\d+|d+",File.readline())
        initialEnergy = float(L[0])
        finalEnergy = float(L[1])
        step = float(L[2])
        X = []
        for i in np.arange(initialEnergy,finalEnergy-0.1,step):
            X.append(i)

        Y = np.loadtxt(txt,skiprows=2)
        Y = savitzky.savitzky_golay(Y,7,3)

        if self.peak_trigger == False:
            if self.toggleShirleyState.get()==1:
                init, end = 0, len(X)

                Ys = shirley(Y,init,end-1)
                y_maxi = Ys.argmax()

                Xs = X-X[y_maxi]
                X = Xs
                Y = Y-min(Y)
                Xs, Ys = Xs[init+1:end], Ys[init+1:end]
                if self.toggleViewState.get() == 1:
                    self.ax.plot(Xs,Ys,'g-')
                if self.toggleGaussianState.get()==1:
                    p = gaussian.fit_iterative(Xs,Ys)
                    A, mu, sigma = p
                    x_n = np.linspace(np.min(Xs),np.max(Xs),100)
                    y = gaussian.gauss(x_n, A=A, mu=mu, sigma=sigma)
                    if self.toggleViewState.get() == 1:
                        self.ax.plot(x_n,y,'r-')
                    if self.run_trigger == True:
                        integrate = np.trapz(y,x_n)
                        txt = txt.split("/")[-1]
                        tp = re.findall(r"[-+]?\d*\.\d+|d+",txt)
                        tp = tp[0].split('.')
                        theta, phi = float(tp[0]), float(tp[1])
                        self.output.write("%f    %f    %f\n"%(integrate,theta,phi))
        else:
            if self.toggleShirleyState.get()==1:
                init, end = 0, len(X)

                Ys = shirley(Y,init,end-1)
                y_maxi = Ys.argmax()

                Xs = X-X[y_maxi]
                X = Xs
                Y = Y-min(Y)
                Xs, Ys = Xs[init+1:end], Ys[init+1:end]

                peak = max(Ys)
                peak_index = np.where(Ys==peak)[0]
                if peak_index>len(Xs)/2:
                    X1, Y1 = Xs[0:peak_index-self.peak_size], Ys[0:peak_index-self.peak_size]
                    X2, Y2 = Xs[peak_index-self.peak_size:-1], Ys[peak_index-self.peak_size:-1]
                else:
                    X1, Y1 = Xs[0:peak_index+self.peak_size], Ys[0:peak_index+self.peak_size]
                    X2, Y2 = Xs[peak_index+self.peak_size:-1], Ys[peak_index+self.peak_size:-1]

                XX = [X1, X2]
                YY = [Y1, Y2]

                if self.toggleViewState.get() == 1:
                    self.ax.plot(Xs,Ys,'g-')
                if self.toggleGaussianState.get()==1:
                    for a in range(len(XX)):
                        p = gaussian.fit_iterative(XX[a],YY[a])
                        A, mu, sigma = p
                        x_n = np.linspace(np.min(Xs),np.max(Xs),100)
                        y = gaussian.gauss(x_n, A=A, mu=mu, sigma=sigma)
                        if self.toggleViewState.get() == 1:
                            self.ax.plot(x_n,y,'r-')
                        if self.run_trigger == True:
                            integrate = np.trapz(y,x_n)
                            txt = txt.split("/")[-1]
                            tp = re.findall(r"[-+]?\d*\.\d+|d+",txt)
                            tp = tp[0].split('.')
                            theta, phi = float(tp[0]), float(tp[1])
                            if a==0:
                                self.output1.write("%f    %f    %f\n"%(integrate,theta,phi))
                            if a==1:
                                self.output2.write("%f    %f    %f\n"%(integrate,theta,phi))
        if self.toggleViewState.get() == 1:
            self.ax.plot(X,Y,'b-')
            self.ax.set_ylim(bottom=0)
            self.canvas.draw()

    def double_peak(self):
        """Activates or deactivates Double Peak"""
        if self.peak_trigger == False:
            self.peak_trigger = True
        else:
            self.peak_trigger = False
        self.plotsimple(self.files[self.index])

    def nextimg(self):
        """Go to the next plot"""
        if self.index < len(self.files)-1:
            self.index += 1
            self.slider.set(self.index)
            self.plotsimple(self.files[self.index])
            self.evlabel["text"] = os.path.basename(self.files[self.index])
            self.label["text"] = ""
        else:
            self.label["text"] = "Reached last image"

    def previmg(self):
        """Go to the previous plot"""
        if self.index == 0:
            self.label["text"] = "Reached first image"
        else:
            self.index -= 1
            self.slider.set(self.index)
            self.plotsimple(self.files[self.index])
            self.evlabel["text"] = os.path.basename(self.files[self.index])
            self.label["text"] = ""

    def increase_psize(self):
        """Increases the peak size in case of double peak."""
        self.peak_size += 1
        self.plotsimple(self.files[self.index])

    def decrease_psize(self):
        """Decreases the peak size in case of double peak."""
        self.peak_size -= 1
        self.plotsimple(self.files[self.index])
 
    def slide_image(self,value):
        """Uses the slider bar to move between plots"""
        if self.Scale:
            self.index = int(value)
            self.plotsimple(self.files[self.index])
            self.evlabel["text"] = os.path.basename(self.files[self.index])
            self.label["text"] = ""

    def config(self):
        """Sets a configuration window."""
        self.config_window = Tkinter.Toplevel(self)
        self.config_frame = Tkinter.Frame(self.config_window)
        self.config_frame.pack(side='top',expand=True, fill='both')
        self.barfig = Figure(figsize=(5,3))
        self.barfig.set_facecolor("#d9d9d9")
        self.ax1 = self.barfig.add_axes([0.1,0,0.1,0.9])
        self.ax2 = self.barfig.add_axes([0.3,0,0.1,0.9])
        self.ax3 = self.barfig.add_axes([0.6,0,0.1,0.9])
        self.ax4 = self.barfig.add_axes([0.8,0,0.1,0.9])
        cb1 = mpl.colorbar.ColorbarBase(self.ax1, cmap = colormaps.afmcold)
        cb2 = mpl.colorbar.ColorbarBase(self.ax2, cmap = colormaps.afmgreen)
        cb3 = mpl.colorbar.ColorbarBase(self.ax3, cmap = colormaps.afmhot)
        cb4 = mpl.colorbar.ColorbarBase(self.ax4, cmap = colormaps.inferno)
        cb1.ax.set_yticklabels([])
        cb2.ax.set_yticklabels([])
        cb3.ax.set_yticklabels([])
        cb4.ax.set_yticklabels([])
        self.barcanvas = FigureCanvasTkAgg(self.barfig,master=self.config_frame)
        self.barcanvas.get_tk_widget().pack(side='top', fill=Tkinter.BOTH, expand=True)
        self.barcanvas.show()
        self.button1 = Tkinter.Button(self.config_frame, text="Blue",command=self.select_blue,width=10)
        self.button2 = Tkinter.Button(self.config_frame, text="Green",command=self.select_green,width=10)
        self.button3 = Tkinter.Button(self.config_frame, text="Red",command=self.select_red,width=10)
        self.button4 = Tkinter.Button(self.config_frame, text="Inferno",command=self.select_inferno,width=10)
        self.button1.pack(side='left')
        self.button2.pack(side='left')
        self.button3.pack(side='left')
        self.button4.pack(side='left')
        self.barcanvas.draw()
        self.checkV_frame = Tkinter.Frame(self.config_window)
        self.checkV_frame.pack(side='bottom',expand=True,fill='both')
        self.checkS_frame = Tkinter.Frame(self.config_window)
        self.checkS_frame.pack(side='bottom',expand=True,fill='both')
        self.checkG_frame = Tkinter.Frame(self.config_window)
        self.checkG_frame.pack(side='bottom',expand=True,fill='both')
        self.toggleView = Tkinter.Checkbutton(self.checkV_frame,variable=self.toggleViewState, onvalue=1,
                                            offvalue=0, text="Toggle graph visualization").pack(side='left')
        self.toggleShirley = Tkinter.Checkbutton(self.checkS_frame,variable=self.toggleShirleyState, onvalue=1,
                                            offvalue=0, text ="Shirley Background Subtraction").pack(side='left')
        self.toggleGaussian = Tkinter.Checkbutton(self.checkG_frame,variable=self.toggleGaussianState, onvalue=1,
                                            offvalue=0, text ="Gaussian Fit").pack(side='left')
        self.pizzaframe = Tkinter.Frame(self.config_window)
        self.pizzaframe.pack(side="top", expand=True, fill="both")
        self.rotatelabel = Tkinter.Label(self.pizzaframe, text="Rotate Pizza (degrees)").pack()
        self.r1 = Tkinter.Entry(self.pizzaframe)
        self.r1.pack(padx=5)
        self.rotatebutton = Tkinter.Button(self.pizzaframe, text="Rotate", command=self.rotate)
        self.rotatebutton.pack()
        self.simmetryframe = Tkinter.Frame(self.config_window)
        self.simmetryframe.pack(side="top", expand=True, fill="both")
        self.simmetrylabel = Tkinter.Label(self.simmetryframe, text="Select system simmetry").pack()
        self.raw = Tkinter.Radiobutton(self.simmetryframe, text = "Raw Data",
                    variable=self.simmetryState, value=0,command=self.change_simmetry).pack()
        self.d60 = Tkinter.Radiobutton(self.simmetryframe, text = u"60\N{DEGREE SIGN}", 
                    variable = self.simmetryState, value=1, command = self.change_simmetry).pack()
        self.d120 = Tkinter.Radiobutton(self.simmetryframe, text = u"120\N{DEGREE SIGN}", 
                    variable = self.simmetryState, value=2, command = self.change_simmetry).pack()
        self.d90 = Tkinter.Radiobutton(self.simmetryframe, text = u"90\N{DEGREE SIGN}",
                    variable = self.simmetryState, value=3, command = self.change_simmetry).pack()
        self.exportframe = Tkinter.Frame(self.config_window)
        self.exportframe.pack(side="top", expand=True, fill="both")
        self.exportbutton = Tkinter.Button(self.exportframe, text="Export Data", command=self.export)
        self.exportbutton.pack()

    def export(self):
        try:
            a = self.files2.split("/")[-1].split(".")
            for tt in range(int(min(self.THETA)),int(max(self.THETA))+1,3):
                for pp in range(int(360*min(self.PHI)/(2*np.pi)),int(360*max(self.PHI)/(2*np.pi))+1,3):
                    p = self.find_nearest(self.gphi*360/(2*np.pi),pp)
                    t = self.find_nearest(self.gtheta,tt)
                    print pp, tt
        except NameError:
            pass

    def find_nearest(self, array, value):
        idx = (np.abs(array-value)).argmin()
        idx2 = np.unravel_index(idx,array.shape)
        return idx2

    def change_simmetry(self):
        if self.simmetryState.get()==0:
            self.plotpizza(rotate=self.degree)
        if self.simmetryState.get()==1:
            self.plotpizza(rotate=self.degree)
        if self.simmetryState.get()==2:
            self.plotpizza(rotate=self.degree)
        if self.simmetryState.get()==3:
            self.plotpizza(rotate=self.degree)

    def select_blue(self,event=None):
        """Selects the afmcold colormap."""
        self.colormap = colormaps.afmcold
        self.plotpizza(rotate = self.degree)

    def select_green(self,event=None):
        """Selects the afmgreen colormap."""
        self.colormap = colormaps.afmgreen
        self.plotpizza(rotate = self.degree)

    def select_red(self,event=None):
        """Selects the afmhot colormap."""
        self.colormap = colormaps.afmhot
        self.plotpizza(rotate = self.degree)

    def select_inferno(self,event=None):
        """Selects the inferno colormap."""
        self.colormap = colormaps.inferno
        self.plotpizza(rotate = self.degree)

    def change_style(self):
        """Set a new window with a Combobox to choose between themes"""
        self.combo_window = Tkinter.Toplevel(self)
        self.combo = ttk.Combobox(self.combo_window, values=self.available_themes)
        self.combo.pack(padx=32, pady=8)
        self.combo.bind('<Return>', self.set_style)
        button = ttk.Button(self.combo_window, text='OK')
        button['command'] = self.set_style
        button.pack(pady=8)

    def set_style(self, event=None):
        """Set the Style to the content of the Combobox"""
        content = self.combo.get()
        try:
            self.style.theme_use(content)
            self.style.configure('TButton', relief='flat')
        except Tkinter.TclError as err:
            tkMessageBox.showerror('Error', err)
        self.combo_window.destroy()

    def run(self):
        """Asks for the name of the output file(s) and go through every plot, 
        integrating the data and saving to the outputs."""
        self.toggleShirleyState.set(1)
        self.toggleGaussianState.set(1)
        self.out_window = Tkinter.Toplevel(self)
        if self.peak_trigger == True:
            Tkinter.Label(self.out_window, text="Output File 1").pack()
            self.e1 = Tkinter.Entry(self.out_window)
            self.e1.pack(padx=5)
            Tkinter.Label(self.out_window, text="Output File2").pack()
            self.e2 = Tkinter.Entry(self.out_window)
            self.e2.pack(padx=5)
        else:
            Tkinter.Label(self.out_window, text="Output File").pack()
            self.e = Tkinter.Entry(self.out_window)
            self.e.pack(padx=5)
        self.okbutton = Tkinter.Button(self.out_window, text="OK",command=self.ok)
        self.okbutton.pack(pady=5)
        self.out_window.bind('<Return>', self.ok)

    def ok(self,event=None):
        """Gets the state of the Label on the Run window and create the output files.
        Creates progress bar to see remaining time."""
        self.run_trigger = True
        self.index = 0
        if self.peak_trigger == True:
            self.output1 = open(self.e1.get()+".xpd",'w')
            self.output1.write("Intensity    Theta    Phi\n")
            self.output2 = open(self.e2.get()+".xpd",'w')
            self.output2.write("Intensity    Theta    Phi\n")
        else:
            self.output = open(self.e.get()+".xpd",'w')
            self.output.write("Intensity    Theta    Phi\n")
        self.out_window.destroy()
        self.progress["value"] = 0
        self.progress["maximum"] = len(self.files)-1
        for i in range(len(self.files)):
            self.plotsimple(self.files[self.index])
            self.evlabel["text"] = os.path.basename(self.files[self.index])
            self.index += 1
            self.progress["value"] = i
            self.update_idletasks()
        if self.peak_trigger==True:
            self.output1.close()
            self.output2.close()
        else:
            self.output.close()
        self.run_trigger = False

    def rotate(self):
        """Rotates Pizza"""
        a = self.r1.get()
        try:
            self.degree = float(a)
            self.plotpizza(rotate=self.degree)
        except ValueError:
            pass


    def makepizza(self):
        """Asks for an output file (searches for a .xpd file, generated in the program.
        The .xpd file must be a file with Intensity, Theta, Phi, in this order.
        Generates the XPD pizza, in this case with 120 degree simmetry, although it can be changed easilly."""
        if self.pizza_trigger==True:
            self.config_window.destroy()
            self.pizza_trigger = False
        self.config()
        self.files2 =  tkFileDialog.askopenfilename(**self.file_opt2)
        self.plotpizza()
        self.pizza_trigger = True

    def plotpizza(self, rotate=0):
        self.ax = self.fig.add_subplot(111)
        self.ax.clear()
        if not self.files2:
            pass
        else:
            self.ext = self.files2.split('.')[-1]
            A = np.loadtxt(self.files2,skiprows=1)

            if self.ext=="xpd":
                intensity, theta, phi = A[:,0], A[:,1], A[:,2]

                intensity = [float(x) for x in intensity]
                theta_degree = theta
                intensity = np.array(intensity)

                theta, phi = np.pi*theta/180., np.pi*phi/180.

                ki = np.array([])
                THETA, PHI  = [], []

                for t in range(int(min(theta_degree)),int(max(theta_degree))+1,3):
                    I_t = intensity[np.where(theta_degree==t)]
                    I_tm = []
                    if self.simmetryState.get()==0:
                        for i in range(len(I_t)):
                            THETA.append(t)
                            PHI.append(3*i)
                        I_m = np.sum(I_t)/len(I_t)
                        NI = (I_t-I_m)/I_m
                        ki = np.append(ki, NI)

                    if self.simmetryState.get()==1:
                        if len(I_t)%4==0:
                            for i in range(len(I_t)/4):
                                I_tm.append((I_t[i]+I_t[i+len(I_t)/4]+I_t[i+2*len(I_t)/4]+I_t[i+3*len(I_t)/4])/2.)
                                THETA.append(t)
                                PHI.append(3*i)
                        else:
                            for i in range(len(I_t)/4+1):
                                if I_t[i]==I_t[-1]:
                                    I_tm.append(I_t[i])
                                    THETA.append(t)
                                    PHI.append(3*i)
                                else:
                                    I_tm.append((I_t[i]+I_t[i+len(I_t)/4]+I_t[i+2*len(I_t)/4]+I_t[i+3*len(I_t)/4])/2.)
                                    THETA.append(t)
                                    PHI.append(3*i)
                        I_m = np.sum(I_tm)/len(I_tm)
                        NI = (I_tm - I_m)/I_m
                        ki = np.append(ki,NI)

                    if self.simmetryState.get()==2:
                        if len(I_t)>230/3.:
                            if len(I_t)%2==0:
                                for i in range(len(I_t)/2):
                                    I_tm.append((I_t[i]+I_t[i+len(I_t)/2])/2.)
                                    THETA.append(t)
                                    PHI.append(3*i)
                            else:
                                for i in range(len(I_t)/2+1):
                                    if I_t[i]==I_t[-1]:
                                        I_tm.append(I_t[i])
                                        THETA.append(t)
                                        PHI.append(3*i)
                                    else:
                                        I_tm.append((I_t[i]+I_t[i+len(I_t)/2])/2.)
                                        THETA.append(t)
                                        PHI.append(3*i)
                        else:
                            for i in range(len(I_t)):
                                I_tm.append(I_t[i])
                                THETA.append(t)
                                PHI.append(3*i)
                        I_m = np.sum(I_tm)/len(I_tm)
                        NI = (I_tm - I_m)/I_m
                        ki = np.append(ki,NI)

                    if self.simmetryState.get()==3:
                        if len(I_t)%3==0:
                            for i in range(len(I_t)/3):
                                I_tm.append((I_t[i]+I_t[i+len(I_t)/3]+I_t[i+2*len(I_t)/3])/2.)
                                THETA.append(t)
                                PHI.append(3*i)
                        else:
                            for i in range(len(I_t)/3+1):
                                if I_t[i]==I_t[-1]:
                                    I_tm.append(I_t[i])
                                    THETA.append(t)
                                    PHI.append(3*i)
                                else:
                                    I_tm.append((I_t[i]+I_t[i+len(I_t)/3]+I_t[i+2*len(I_t)/3])/2.)
                                    THETA.append(t)
                                    PHI.append(3*i)
                        I_m = np.sum(I_tm)/len(I_tm)
                        NI = (I_tm - I_m)/I_m
                        ki = np.append(ki,NI)

            if self.ext=="hol":
                THETA, PHI, ki = A[:,0], A[:,1], A[:,4]
                THETA = np.array(THETA)
                ki = np.array(ki)

            if self.ext=="dat":
                THETA, PHI, ki = A[:,0], A[:,1], A[:,2]
                THETA = np.array(THETA)
                ki = np.array(ki)

            PHI = np.array(PHI)
            #Shift to get the XPD pizza centralized
            PHI = np.pi*PHI/180 - np.pi*rotate/180

            self.PHI_ = PHI
            self.THETA_ = THETA
            self.ki_ = ki

            if self.simmetryState.get()==1:
                for i in range(1,6):
                    PHI = np.append(PHI, self.PHI_ + i*60*np.pi/180)
                    THETA = np.append(THETA,self.THETA_)
                    ki = np.append(ki, self.ki_)
            if self.simmetryState.get()==2:
                for i in range(1,3):
                    PHI = np.append(PHI, self.PHI_ + i*120*np.pi/180)
                    THETA = np.append(THETA, self.THETA_)
                    ki = np.append(ki, self.ki_)
            if self.simmetryState.get()==3:
                for i in range(1,4):
                    PHI = np.append(PHI, self.PHI_ + i*90*np.pi/180)
                    THETA = np.append(THETA, self.THETA_)
                    ki = np.append(ki, self.ki_)

            PHI = PHI - min(PHI)  #Shifts PHI to zero
            del_index = np.where(PHI>2*np.pi)  #Index of everything greater than 2pi
            THETA = np.delete(THETA, del_index)  #Deletes indexes
            ki = np.delete(ki, del_index)
            PHI = np.delete(PHI, del_index)

            self.THETA, self.PHI = THETA, PHI

            if self.simmetryState.get()==0:
                gx, gy = np.mgrid[min(THETA):max(THETA):600j, min(PHI):max(PHI):600j]
                self.z = interpolate.griddata((THETA,PHI),ki,(gx,gy), method='cubic',fill_value=ki.min())
                self.gphi, self.gtheta = gy, gx
            else:
                xx = THETA*np.cos(PHI)
                yy = THETA*np.sin(PHI)

                gx, gy = np.mgrid[min(xx):max(xx):600j, min(yy):max(yy):600j]

                self.z = interpolate.griddata((xx, yy),ki,(gx,gy), method='cubic',fill_value=ki.min())
                self.gphi = np.arctan2(gy, gx)
                self.gtheta = np.sqrt(gx**2 + gy**2)
 
                self.gphi = self.gphi -np.amin(self.gphi)

                h, l = np.where(self.gtheta<min(THETA))
                self.z[h,l] = min(ki)

            norm = mpl.colors.Normalize(self.z.flatten().min(),self.z.flatten().max())

            self.fig.delaxes(self.ax)
            self.ax2 = self.fig.add_subplot(1,1,1, polar = True)
            self.ax2.pcolormesh(self.gphi, self.gtheta, self.z,cmap=self.colormap,norm=norm)
            self.ax2.set_yticks([min(THETA),max(THETA)])
            self.ax2.grid(True)
            self.ax2.set_rmax(max(THETA))
            self.fig.savefig(self.files2+".jpg")
            self.canvas.draw()
            self.fig.delaxes(self.ax2)

MainWindow().mainloop()
