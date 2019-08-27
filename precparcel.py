#/usr/env/python
import matplotlib
matplotlib.use('wx')
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons
from parcelmod import *
import time
from colorsys import *

rcParams['mathtext.default']='regular'
ax = subplot(111)
labels=['s','mse','ths','thl','mse app','thl app','thl lin/ice','thq palch']

tinit=20.0
texcinit=0.0
qinit=0.0
rhinit=85
lapseinit=6.5
epsinit=0.0
ipiso=0
cmse=1
iice=1
isatice=1
ismooth=1
ipz=0
iprecliq=0
iprecice=0
plinit=1500
piinit=1500

subplots_adjust(left=0.25, bottom=0.42)
axcolor = 'lightgoldenrodyellow'

axt=[]
for i in range(8):
    axt.append(axes([0.25, 0.16+0.025*i, 0.65, 0.02]))

tsl=[]
tsl.append(Slider(axt[0], 'Tv  surf (C)', -20.0, 50.0, valinit=tinit))
tsl.append(Slider(axt[1], 'T  surf exc (C)', -2.0, 5.0, valinit=texcinit))
tsl.append(Slider(axt[2], 'qt exc (g/kg)', -2.0, 4.0, valinit=qinit))
tsl.append(Slider(axt[3], 'Env RH (%)', 0.0, 100.0, valinit=rhinit))
tsl.append(Slider(axt[4], 'Trop T_v lapse', 0.0, 10.0, valinit=lapseinit))
tsl.append(Slider(axt[5], 'Ent constant *1e-3', 0.0, 2.0, valinit=epsinit))
tsl.append(Slider(axt[6], 'Liq prec length', 0.0, 5000.0, valinit=plinit))
tsl.append(Slider(axt[7], 'Ice prec length', 0.0, 5000.0, valinit=piinit))

def update(val):
    tvsurf=tsl[0].val
    tvsurfexc=tsl[1].val
    qtsurfexc=tsl[2].val
    envrh=tsl[3].val
    lapse=tsl[4].val
    eps  =tsl[5].val
    plenl =tsl[6].val
    pleni =tsl[7].val
    precparcel.entrain(tvsurf+273.15,tvsurfexc,qtsurfexc/1000.,envrh/100.,-lapse/1000.,eps,ipiso,cmse,iice,isatice,ismooth,ipz,iprecliq,iprecice,plenl,pleni)
    ax.clear()
    # plot stuff
    ax.set_color_cycle(['#CC0000','#FF6600','#FFCC00','#99FF33','#66FFFF','#0000FF','#000000','#663300','#9900CC'])
    ax.plot(parceldata.tvent-parceldata.tvenv,parceldata.z,parceldata.tvmse-parceldata.tvenv,parceldata.z,parceldata.tvths-parceldata.tvenv,parceldata.z,parceldata.tvthl-parceldata.tvenv,parceldata.z,parceldata.tvmseapr-parceldata.tvenv,parceldata.z,parceldata.tvthlapr-parceldata.tvenv,parceldata.z,parceldata.tvthllin-parceldata.tvenv,parceldata.z,parceldata.tvtheapr-parceldata.tvenv,parceldata.z)
    #ax.plot(parceldata.tvent-parceldata.tvenv,parceldata.z,parceldata.tvmse-parceldata.tvenv,parceldata.z)
    ax.legend(labels,ncol=3)
    ax.set_title(r'Delta Tv')
    ax.set_xlim(-10,10)
    ax.set_ylim(0,20000)
    draw()

for i in range(8):
    tsl[i].on_changed(update)

resetax = axes([0.8, 0.01, 0.1, 0.07])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    for i in range(8):
        tsl[i].reset()

button.on_clicked(reset)

def pisoupdate(val):
    global ipiso
    ipiso=(ipiso+1)%2
    update(val)

def cmseupdate(val):
    global cmse
    cmse=(cmse+1)%2
    update(val)

def isaticeupdate(val):
    global isatice
    isatice=(isatice+1)%2
    update(val)

def iiceupdate(val):
    global iice
    iice=(iice+1)%2
    update(val)

def ismoothupdate(val):
    global ismooth
    ismooth=(ismooth+1)%2
    update(val)

def ipzupdate(val):
    global ipz
    ipz=(ipz+1)%2
    update(val)

def iprecliqupdate(val):
    global iprecliq
    iprecliq=(iprecliq+1)%2
    update(val)
                
def ipreciceupdate(val):
    global iprecice
    iprecice=(iprecice+1)%2
    update(val)

               
#pp = axes([0.6, 0.01, 0.2, 0.05])
#check = CheckButtons(pp, ['p-iso'], [False])
#check.on_clicked(pisoupdate)

pp = axes([0.6, 0.01, 0.2, 0.07])
check = CheckButtons(pp, ['Entrain ~1/z'], [False])
check.on_clicked(ipzupdate)

pp2 = axes([0.4, 0.01, 0.2, 0.07])
check2 = CheckButtons(pp2, ['MSE correction'], [True])
check2.on_clicked(cmseupdate)

pp3 = axes([0.2, 0.01, 0.2, 0.07])
check3 = CheckButtons(pp3, ['Ice'], [True])
check3.on_clicked(iiceupdate)

pp4 = axes([0.0, 0.01, 0.2, 0.07])
check4 = CheckButtons(pp4, ['Smooth freeze'], [True])
check4.on_clicked(ismoothupdate)

pp5 = axes([0.0, 0.08, 0.2, 0.07])
check5 = CheckButtons(pp5, ['Prec liq'], [False])
check5.on_clicked(iprecliqupdate)

pp6 = axes([0.2, 0.08, 0.2, 0.07])
check6 = CheckButtons(pp6, ['Prec ice'], [False])
check6.on_clicked(ipreciceupdate)

pp7 = axes([0.4, 0.08, 0.2, 0.07])
check7 = CheckButtons(pp7, ['Satur ice'], [True])
check7.on_clicked(isaticeupdate)

time.sleep(0.2)
show()
