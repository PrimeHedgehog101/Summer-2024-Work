############################ Packages ###############################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from Stream_arrows import streamQuiver

############################ Variables #################################

mu_0 = 4*np.pi*10
mu_d = 2*mu_0
mu_b = mu_0
eta_1 = 80
eta_2 = 50
H_inf = 1
R_0 = 1
a = 1
b = 2
N = 1000

############################# Global set-up ################################

x = np.linspace(-2,2,N)
y = np.linspace(-2,2,N)
X , Y = np.meshgrid(x,y)

R = np.sqrt(X**2+Y**2)

############################ Functions ###############################

def ferro1():
    
    A = (2*mu_0 * H_inf)/(mu_0+mu_d)
    B = H_inf
    D = (mu_0-mu_d)/(mu_0+mu_d) * R_0**2 * H_inf
    
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    
    indrop = R < R_0
    outdrop = R >= R_0 
    
    x1 = X[indrop]
    x2 = X[outdrop]
    y1 = Y[indrop]
    y2 = Y[outdrop]
    
    U[indrop] = 0 # i component for r<R_0, (X[indrop] or Y[indrop] only)
    V[indrop] = -A # j component for r<R_0
    
    U[outdrop] = 2*D*x2*y2/(x2**2+y2**2)**2 # i component for r>R_0
    V[outdrop] = -B-(D/(x2**2+y2**2))+(2*D*y2**2/(x2**2+y2**2)**2) # j component for r>R_0
    
    U = np.nan_to_num(U)
    V = np.nan_to_num(V)
    
    return U,V


def plot1():
    
    U,V = ferro1()[0] , ferro1()[1]
    
    fig = plt.figure(figsize=(7,6),dpi=200)
    gs = GridSpec(1,3,width_ratios=[1,1,0.2],wspace=None,hspace=None)
    
    
    ax1 = fig.add_subplot(gs[0,0])
    
    x_points = np.linspace(x[0],0,20) #20 with minlength=0.8, density=5
    y_points = np.linspace(y[N-1],y[N-1],20)
    
    stream_points = np.column_stack((x_points,y_points))
    ax1.streamplot(X,Y,U,V,start_points=stream_points,color='black',arrowstyle='->',linewidth=1,minlength=0,density=5)
    ax1.set_xlim(x[0],0)
    ax1.set_ylim(y[0],y[N-1])
    ax1.set_aspect('auto')
    circle1 = plt.Circle((0,0),R_0,color='black',fill=False,linewidth=2)
    ax1.add_patch(circle1)
    
    
    ax2 = fig.add_subplot(gs[0,1])
    Magnitude = np.sqrt(U**2+V**2)
    contour = ax2.contourf(X,Y,Magnitude,cmap='viridis',levels=15)
    ax2.set_xlim(0,x[N-1])
    ax2.set_ylim(y[0],y[N-1])
    ax2.set_aspect('auto')
    ax2.yaxis.set_visible(False)
    
    
    circle2 = plt.Circle((0,0),R_0,color='black',fill=False,linewidth=2)
    ax2.add_patch(circle2)
    
    ax3 = fig.add_subplot(gs[0,2])
    fig.colorbar(contour,cax=ax3,orientation='vertical')
    ax3.set_aspect(20)
    
    
    plt.subplots_adjust(left=0,right=1,top=1,bottom=0,wspace=0,hspace=0)
    
    
    plt.show()
    #fig.savefig('ferro_1a.pdf',bbox_inches='tight')



def ferro2():
    
    alpha = (b**2-a**2)*mu_0**2+(a**2+b**2)*mu_0*(mu_d+mu_b)+(b**2-a**2)*mu_d*mu_b
    A = 4*mu_0*mu_b*H_inf*b**2 / alpha
    B = 2*mu_b*(mu_d+mu_0)*H_inf*b**2 / alpha
    C = 2*mu_b*(mu_0-mu_d)*H_inf*b**2*a**2 / alpha
    D = ((a**2-b**2)*mu_0**2+(a**2+b**2)*mu_0*(mu_b-mu_d)+(b**2-a**2)*mu_d*mu_b)*H_inf*b**2 / alpha
    
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    
    indrop = R < a
    outdrop = (R>a) & (R<b)
    wall = R > b
    
    x1 = X[indrop]
    x2 = X[outdrop]
    x3 = X[wall]
    y1 = Y[indrop]
    y2 = Y[outdrop]
    y3 = Y[wall]
    
    U[indrop] = 0
    V[indrop] = -A
    
    U[outdrop] = 2*C*x2*y2/(x2**2+y2**2)**2
    V[outdrop] = -B-(C/(x2**2+y2**2))+(2*C*y2**2/(x2**2+y2**2)**2)
    
    U[wall] = 2*D*x3*y3/(x3**2+y3**2)**2
    V[wall] = -H_inf-(D/(x3**2+y3**2))+(2*D*y3**2/(x3**2+y3**2)**2)
    
    U = np.nan_to_num(U)
    V = np.nan_to_num(V)
    
    return U,V

def plot2():
    
    U,V = ferro2()[0] , ferro2()[1]
    
    fig = plt.figure(figsize=(7,6),dpi=200)
    gs = GridSpec(1,3,width_ratios=[1,1,0.2],wspace=None,hspace=None)
    
    
    ax1 = fig.add_subplot(gs[0,0])
    
    x_points = np.linspace(x[0],0,20) #20 with minlength=0.8, density=5
    y_points = np.linspace(y[N-1],y[N-1],20)
    
    stream_points = np.column_stack((x_points,y_points))
    ax1.streamplot(X,Y,U,V,start_points=stream_points,color='black',arrowstyle='->',linewidth=1,minlength=0,density=5)
    ax1.set_xlim(x[0],0)
    ax1.set_ylim(y[0],y[N-1])
    ax1.set_aspect('auto')
    
    circle1a = plt.Circle((0,0),a,color='black',fill=False,linewidth=2)
    ax1.add_patch(circle1a)
    
    circle1b = plt.Circle((0,0),b,color='black',fill=False,linewidth=2)
    ax1.add_patch(circle1b)
    
    
    ax2 = fig.add_subplot(gs[0,1])
    Magnitude = np.sqrt(U**2+V**2)
    contour = ax2.contourf(X,Y,Magnitude,cmap='viridis',levels=15)
    ax2.set_xlim(0,x[N-1])
    ax2.set_ylim(y[0],y[N-1])
    ax2.set_aspect('auto')
    ax2.yaxis.set_visible(False)
    
    
    circle2a = plt.Circle((0,0),a,color='black',fill=False,linewidth=2)
    ax2.add_patch(circle2a)
    
    circle2b = plt.Circle((0,0),b,color='black',fill=False,linewidth=2)
    ax2.add_patch(circle2b)
    
    ax3 = fig.add_subplot(gs[0,2])
    fig.colorbar(contour,cax=ax3,orientation='vertical')
    ax3.set_aspect(20)
    
    
    plt.subplots_adjust(left=0,right=1,top=1,bottom=0,wspace=0,hspace=0)
    
    
    plt.show()
    fig.savefig('ferro_6a.pdf',bbox_inches='tight')


def ferro3():
    
    const = (H_inf**2 * mu_0 * (mu_0-mu_d)**2) / (8*(mu_0+mu_d)**2*eta_2*(eta_1-2*eta_2))
    
    U_r = np.zeros_like(X)
    V_t = np.zeros_like(Y)
    U_x = np.zeros_like(X)
    V_y = np.zeros_like(Y)
    
    indrop = R < R_0
    outdrop = R >= R_0 
   
    x1 = X[indrop]
    x2 = X[outdrop]
    y1 = Y[indrop]
    y2 = Y[outdrop]
    r_1 = np.sqrt(X[indrop]**2+Y[indrop]**2)
    r_2 = np.sqrt(X[outdrop]**2+Y[outdrop]**2)
    
    U_r[indrop] = (-3*(eta_1-3*eta_2)*r_1+((eta_1-5*eta_2)*r_1**3/a**2))*const*(x1**2-y1**2)/(x1**2+y1**2) # i component for r<R_0, (X[indrop] or Y[indrop] only)
    V_t[indrop] = (-3*(eta_1-3*eta_2)*r_1+(2*(eta_1-5*eta_2)*r_1**3/a**2))*-1*const*(2*x1*y1)/(x1**2+y1**2) # j component for r<R_0
    
    U_r[outdrop] = (-3*(eta_1-eta_2)*r_2**(-1)+a**2*(eta_1+eta_2)*r_2**(-3))*const*a**2*(x2**2-y2**2)/(x2**2+y2**2)  # i component for r>R_0
    V_t[outdrop] = (eta_1+eta_2)*r_2**(-3)*const*a**4*(2*x2*y2)/(x2**2+y2**2) # j component for r>R_0
    
    U_r = np.nan_to_num(U_r)
    V_t = np.nan_to_num(V_t)
    
    U_x[indrop] = (x1*U_r[indrop]-y1*V_t[indrop])/r_1
    U_x[outdrop] = (x2*U_r[outdrop]-y2*V_t[outdrop])/r_2
    
    V_y[indrop] = (r_1*V_t[indrop]+y1*U_x[indrop])/x1
    V_y[outdrop] = (r_2*V_t[outdrop]+y2*U_x[outdrop])/x2

    return U_x , V_y 


def plot3():
    
    U,V = ferro3()[0] , ferro3()[1]
    
    fig = plt.figure(figsize=(7,6),dpi=200)
    ax=fig.add_subplot(111)
    #gs = GridSpec(1,3,width_ratios=[1,1,0.2],wspace=None,hspace=None)
    
    
    #ax1 = fig.add_subplot(gs[0,0])
    
    x_points = np.linspace(x[0],x[-1],20) #20 with minlength=0.8, density=5
    y_points = np.linspace(y[N-1],y[N-1],20)
    
    theta_c=np.linspace(0,2*np.pi,100)
    xc=a*np.cos(theta_c)
    yc=a*np.sin(theta_c)
    
    theta_s=np.linspace(0,2*np.pi,41)+np.pi/41
    #theta_s=theta_s[:-1]
    xs=a*np.cos(theta_s)
    ys=a*np.sin(theta_s)
    points_s = list(zip(xs,ys))
    
    #x_values = [i[0] for i in points_s]
    #y_values = [i[1] for i in points_s]
    
    #plt.scatter(x_values,y_values,marker='o')
    
    stream_points = points_s
    stream = plt.streamplot(X,Y,U,V,start_points=stream_points,color='black',arrowstyle='-',linewidth=1,minlength=0,maxlength=0.7,density=5,broken_streamlines=False)
    
    plt.plot(xc,yc,color='black',linewidth=2)
    
    
    qv = streamQuiver(ax,stream,spacing=0.8, scale=40)

    
    plt.show()
    #fig.savefig('ferro_5a.pdf',bbox_inches='tight')


def plot3b():
    
    U,V = ferro3()[0] , ferro3()[1]
    
    fig = plt.figure(figsize=(7,6),dpi=200)
    #ax=fig.add_subplot(111)
    gs = GridSpec(1,3,width_ratios=[1,1,0.2],wspace=None,hspace=None)
    
    
    ax1 = fig.add_subplot(gs[0,0])
    
    #x_points = np.linspace(x[0],x[-1],20) #20 with minlength=0.8, density=5
    #y_points = np.linspace(y[N-1],y[N-1],20)
    
    
    
    theta_s=np.linspace(0,2*np.pi,41)+np.pi/41
    #theta_s=theta_s[:-1]
    xs=a*np.cos(theta_s)
    ys=a*np.sin(theta_s)
    points_s = list(zip(xs,ys))
    
    x_points = np.linspace(x[0],-0.05,14) #20 with minlength=0.8, density=5
    y_points = np.linspace(y[N-1],y[N-1],14)
    
    stream_points = np.column_stack((x_points,y_points))
    stream2 = ax1.streamplot(X,Y,U,V,start_points=stream_points,color='black',arrowstyle='-',linewidth=1,minlength=0,maxlength=5,density=5,broken_streamlines=False)
    qv = streamQuiver(ax1,stream2,spacing=0.2, scale=15)
    
    x_points2 = np.linspace(x[0],-0.05,14) #20 with minlength=0.8, density=5
    y_points2 = np.linspace(y[0],y[0],14)
    
    stream_points2 = np.column_stack((x_points2,y_points2))
    stream3 = ax1.streamplot(X,Y,U,V,start_points=stream_points2,color='black',arrowstyle='-',linewidth=1,minlength=0,maxlength=5,density=5,broken_streamlines=False)
    qv = streamQuiver(ax1,stream3,spacing=0.2, scale=15)
    #x_values = [i[0] for i in points_s]
    #y_values = [i[1] for i in points_s]
    
    #plt.scatter(x_values,y_values,marker='o')
    
    stream_points = points_s
    
    ax1.streamplot(X,Y,U,V,start_points=stream_points,color='black',arrowstyle='-',linewidth=1,minlength=0,maxlength=1,density=5,broken_streamlines=False)
    ax1.set_xlim(x[0],0)
    ax1.set_ylim(y[0],y[N-1])
    ax1.set_aspect('auto')
    circle1 = plt.Circle((0,0),R_0,color='black',fill=False,linewidth=2)
    ax1.add_patch(circle1)
    
    ax2 = fig.add_subplot(gs[0,1])
    Magnitude = np.sqrt(U**2+V**2)
    contour = ax2.contourf(X,Y,Magnitude,cmap='viridis',levels=15)
    ax2.set_xlim(0,x[N-1])
    ax2.set_ylim(y[0],y[N-1])
    ax2.set_aspect('auto')
    ax2.yaxis.set_visible(False)
    
    circle2 = plt.Circle((0,0),R_0,color='black',fill=False,linewidth=2)
    ax2.add_patch(circle2)
    
    ax3 = fig.add_subplot(gs[0,2])
    fig.colorbar(contour,cax=ax3,orientation='vertical')
    ax3.set_aspect(20)
    
    plt.subplots_adjust(left=0,right=1,top=1,bottom=0,wspace=0,hspace=0)
    

    stream = ax1.streamplot(X,Y,U,V,start_points=stream_points,color='black',arrowstyle='-',linewidth=1,minlength=0,maxlength=0.24,density=5,broken_streamlines=False)
    qv = streamQuiver(ax1,stream,spacing=0.2, scale=15)

    
    plt.show()
    #fig.savefig('ferro_5b.pdf',bbox_inches='tight')









#################################### Back- up ##################################

#circle = plt.Circle((0,0),R_0,color='black',fill=False,linewidth=1)
# plt.figure(dpi=200)
# plt.quiver(X,Y,U,V)
# plt.gca().add_artist(circle)
# plt.xlim([-2,2])
# plt.ylim([-2,2])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

