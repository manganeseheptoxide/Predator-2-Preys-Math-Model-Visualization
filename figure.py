from scipy.integrate import odeint
from scipy.integrate import OdeSolution
from math import log as log
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation

# params from the paper

# i put it outside the ODE function so i can easily find it when i need to edit,
# can prolly remove about 200 lines of code if i did not redeclare it inside the functions

r_0 = log(3)
m_0 = 2/7
cVV_0 = log(3)-(2/7)
e_0 = 0.91

s_0 = log(4.5)
n_0 = 4/5
cSS_0 = (log(4.5)/100)-(4/500)

u_0 = log(5)
p_0 = 2/11
cLL_0 = (log(5)/30)-(2/330)

VSL_0 = [0.1, 2, 1]
a = 0.2
b = 3
years = 10

def dtChange(VSL_vector, t, a, b):
    #fox param
    r = r_0
    m = m_0
    cVV = cVV_0
    e = e_0
    #cottontail param
    s = s_0
    n = n_0
    cSS = cSS_0
    #hare param
    u = u_0
    p = p_0
    cLL = cLL_0

    V = VSL_vector[0]
    S = VSL_vector[1]
    L = VSL_vector[2]

    dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
    dSdt = S*(s - cSS*S - n - a*V)
    dLdt = L*(u - cLL*L - p - b*V)

    return [dVdt, dSdt, dLdt]

t = np.linspace(0, years, 1000)
x = odeint(dtChange, VSL_0, t, args=(a, b))

V_data = x[:,0]
S_data = x[:,1]
L_data = x[:,2]

#equlibrium stuff

#eq_fig = plt.figure()
#fig_ax = eq_fig.subplots()
#plt.subplots_adjust(top = 10, bottom = 9)

###feasib

#E0_f = 'N/A'
#E1_f = 'Feasible' if m_0 <= r_0 else 'Not Feasible'
# E2_f = 'Feasible' if n_0 <= s_0 else 'Not Feasible'
# E3_f = 'Feasible' if a*e_0*s_0 + cSS_0*r_0 >= a*e_0*n_0 + cSS_0*m_0 and a*m_0 + cVV_0*s_0 >= a*r_0 + cVV_0*n_0 else 'Not Feasible'
# E4_f = 'Feasible' if p_0 <= u_0 else 'Not Feasible'
#E5_f = 'Feasible' if b*e_0*u_0 + cLL_0*r_0 >= b*e_0*p_0 + cLL_0*m_0 and b*m_0 + cVV_0*u_0 >= b*r_0 + cVV_0*p_0 else 'Not Feasible'
#E6_f = 'Feasible' if n_0 <= s_0 and p_0 <= u_0 else 'Not Feasible'

#con1 = b*e_0*cSS_0*u_0 + a*e_0*cLL_0*s_0 + cSS_0*cLL_0*r_0 
#con2 = b*e_0*cSS_0*p_0 + a*e_0*cLL_0*n_0 + cSS_0*cLL_0*m_0*a*b*e_0*p_0 + (b**2)*e_0*s_0 + a*cLL_0*m_0 + cVV_0*cLL_0*s_0
#con3 = a*b*e_0*u_0 + (b**2)*e_0*n_0 + a*cLL_0*r_0 + cVV_0*cLL_0*n_0*(a**2)*e_0*u_0 + a*b*e_0*n_0 + b*cSS_0*m_0 + cVV_0*cSS_0*u_0
#con4 = (a**2)*e_0*p_0 + a*b*e_0*s_0 + b*cSS_0*r_0 + cVV_0*cSS_0*p_0 

#E7_f = 'Feasible' if con1 >= con2 and con2 >= con3 and con3 >= con4 else 'Not Feasible'

#table_content = [[E0_f, ""],
#                 [E1_f, ""],
#                 [E2_f, ""],
#                 [E3_f, ""],
#                 [E4_f, ""],
#                 [E5_f, ""],
#                 [E6_f, ""],
#                 [E7_f, ""]]

#row_labels = ['E0', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7']
#col_labels = ['Equilibrium', 'Feasibility', 'Stability']

#table = fig_ax.table(cellText=table_content,
#                  colWidths = [0.4]*3,
#                  rowLabels = row_labels,
#                  colLabels = col_labels,
#                  loc='center')


#back parameter sliders

back_param_slider = plt.figure()
param = back_param_slider.subplots()
back_param_slider.clear()

r_ = back_param_slider.add_axes([0.2, 0.95, 0.65, 0.03])
m_ = back_param_slider.add_axes([0.2, 0.9, 0.65, 0.03])
cVV_ = back_param_slider.add_axes([0.2, 0.85, 0.65, 0.03])
e_ = back_param_slider.add_axes([0.2, 0.8, 0.65, 0.03])

slider_r = Slider(r_, "r", valmin = 0.0, valmax = 2*r_0, valinit = r_0, valstep = 0.0001)
slider_m = Slider(m_, "m", valmin = 0.0, valmax = 2*m_0, valinit = m_0, valstep = 0.0001)
slider_cVV = Slider(cVV_, "cVV", valmin = 0.0, valmax = 4*cVV_0, valinit = cVV_0, valstep = 0.00001)
slider_e = Slider(e_, "e", valmin = 0.0, valmax = 2*e_0, valinit = e_0, valstep = 0.001)

s_ = back_param_slider.add_axes([0.2, 0.7, 0.65, 0.03])
n_ = back_param_slider.add_axes([0.2, 0.65, 0.65, 0.03])
cSS_ = back_param_slider.add_axes([0.2, 0.6, 0.65, 0.03])

slider_s = Slider(s_, "s", valmin = 0.0, valmax = 2*s_0, valinit = s_0, valstep = 0.0001)
slider_n = Slider(n_, "n", valmin = 0.0, valmax = 2*n_0, valinit = n_0, valstep = 0.0001)
slider_cSS = Slider(cSS_, "cSS", valmin = 0.0, valmax = 4*cSS_0, valinit = cSS_0, valstep = 0.00001)

u_ = back_param_slider.add_axes([0.2, 0.5, 0.65, 0.03])
p_ = back_param_slider.add_axes([0.2, 0.45, 0.65, 0.03])
cLL_ = back_param_slider.add_axes([0.2, 0.4, 0.65, 0.03])

slider_u = Slider(u_, "u", valmin = 0.0, valmax = 2*u_0, valinit = u_0, valstep = 0.0001)
slider_p = Slider(p_, "p", valmin = 0.0, valmax = 2*p_0, valinit = p_0, valstep = 0.0001)
slider_cLL = Slider(cLL_, "cLL", valmin = 0.0, valmax = 4*cLL_0, valinit = cLL_0, valstep = 0.00001)

reset_b1 = back_param_slider.add_axes([0.4, 0.325, 0.3, 0.04])
reset_back_param = Button(reset_b1, 'Reset Back Parameters', color='gold',hovercolor='skyblue')

V_= back_param_slider.add_axes([0.2, 0.25, 0.65, 0.03])
S_= back_param_slider.add_axes([0.2, 0.2, 0.65, 0.03])
L_= back_param_slider.add_axes([0.2, 0.15, 0.65, 0.03])

slider_V = Slider(V_, "V", valmin = 0.0, valmax = 20, valinit = VSL_0[0], valstep = 0.1)
slider_S = Slider(S_, "S", valmin = 0.0, valmax = 20, valinit = VSL_0[1], valstep = 0.1)
slider_L = Slider(L_, "L", valmin = 0.0, valmax = 20, valinit = VSL_0[2], valstep = 0.1)

E_0 = back_param_slider.add_axes([0.05, 0.9, 0.05, 0.04])
E0_button = Button(E_0, 'E0', color='gold',hovercolor='skyblue')
E_1 = back_param_slider.add_axes([0.05, 0.85, 0.05, 0.04])
E1_button = Button(E_1, 'E1', color='gold',hovercolor='skyblue')
E_2 = back_param_slider.add_axes([0.05, 0.8, 0.05, 0.04])
E2_button = Button(E_2, 'E2', color='gold',hovercolor='skyblue')
E_3 = back_param_slider.add_axes([0.05, 0.75, 0.05, 0.04])
E3_button = Button(E_3, 'E3', color='gold',hovercolor='skyblue')
E_4 = back_param_slider.add_axes([0.05, 0.7, 0.05, 0.04])
E4_button = Button(E_4, 'E4', color='gold',hovercolor='skyblue')
E_5 = back_param_slider.add_axes([0.05, 0.65, 0.05, 0.04])
E5_button = Button(E_5, 'E5', color='gold',hovercolor='skyblue')
E_6 = back_param_slider.add_axes([0.05, 0.6, 0.05, 0.04])
E6_button = Button(E_6, 'E6', color='gold',hovercolor='skyblue')
E_7 = back_param_slider.add_axes([0.05, 0.55, 0.05, 0.04])
E7_button = Button(E_7, 'E7', color='gold',hovercolor='skyblue')

reset_b3 = back_param_slider.add_axes([0.4, 0.075, 0.3, 0.04])
reset_population = Button(reset_b3, 'Reset Population', color='gold',hovercolor='skyblue')

#simulation figure

sim = plt.figure()
ax = sim.subplots()
plt.subplots_adjust(bottom = 0.4)
ax.axis([0, 10, 0, 10]) #########

V_plt, = ax.plot(t, V_data, 'r--')
S_plt, = ax.plot(t, S_data, 'blue')
L_plt, = ax.plot(t, L_data, 'black')
ax.legend(["Fox", "Cottontail", "Hare"])
ax.set_ylabel('Population')
ax.set_xlabel('time')

ax_slide_a = sim.add_axes([0.25, 0.25, 0.65, 0.03])
ax_slide_b = sim.add_axes([0.25, 0.20, 0.65, 0.03])
ax_slide_anb = sim.add_axes([0.25, 0.15, 0.65, 0.03])

slider_a = Slider(ax_slide_a, "a", valmin = 0.0, valmax = 6 if a < 6 else a*2, valinit = a, valstep = 0.01)
slider_b = Slider(ax_slide_b, "b", valmin = 0.0, valmax = 6 if b < 6 else b*2, valinit = b, valstep = 0.01)
slider_anb = Slider(ax_slide_anb, "a & b", valmin = 0.0, valmax = 6 if b < 6 else b*2, valinit = b, valstep = 0.01)

reset_b2 = sim.add_axes([0.5, 0.05, 0.3, 0.04])
reset_bifurcation = Button(reset_b2, 'Reset Parameters', color='gold',hovercolor='skyblue')



def update(val):

    #bifurcation parameters
    a_cur = slider_a.val
    b_cur = slider_b.val

    #back parameters
    r_0 = slider_r.val
    m_0 = slider_m.val
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = slider_s.val
    n_0 = slider_n.val
    cSS_0 = slider_cSS.val

    u_0 = slider_u.val
    p_0 = slider_p.val
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

def E0(val): 

    #bifurcation parameters
    a_cur = slider_a.val
    b_cur = slider_b.val

    #back parameters
    r_0 = 0.1
    m_0 = 0.2
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = 0.1
    n_0 = 0.2
    cSS_0 = slider_cSS.val

    u_0 = 0.1
    p_0 = 0.2
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

E0_button.on_clicked(E0)

def E1(val): 

    #bifurcation parameters
    a_cur = 1
    b_cur = 3

    #back parameters
    r_0 = slider_r.val
    m_0 = slider_m.val
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = slider_s.val
    n_0 = slider_n.val
    cSS_0 = slider_cSS.val

    u_0 = slider_u.val
    p_0 = slider_p.val
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

E1_button.on_clicked(E1)

def E2(val): 

    #bifurcation parameters
    a_cur = 0.1
    b_cur = slider_b.val

    #back parameters
    r_0 = 0
    m_0 = 1
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = 0.2
    n_0 = 0.1
    cSS_0 = slider_cSS.val

    u_0 = 0.1
    p_0 = 0.2
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

E2_button.on_clicked(E2)

def E4(val):

    #bifurcation parameters
    a_cur = 1
    b_cur = 0

    #back parameters
    r_0 = 0.1
    m_0 = 0.2
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = 0.1
    n_0 = 0.2
    cSS_0 = slider_cSS.val

    u_0 = 0.2
    p_0 = 0.1
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

E4_button.on_clicked(E4)

def E3(val):

    #bifurcation parameters
    a_cur = slider_a.val
    b_cur = 5

    #back parameters
    r_0 = slider_r.val
    m_0 = slider_m.val
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = slider_s.val
    n_0 = slider_n.val
    cSS_0 = slider_cSS.val

    u_0 = slider_u.val
    p_0 = slider_p.val
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

E3_button.on_clicked(E3)

def E5(val):

    #bifurcation parameters
    a_cur = slider_a.val
    b_cur = 0.4

    #back parameters
    r_0 = slider_r.val
    m_0 = slider_m.val
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = 0.1
    n_0 = 0.2
    cSS_0 = slider_cSS.val

    u_0 = slider_u.val
    p_0 = slider_p.val
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

E5_button.on_clicked(E5)

def E6(val):

    #bifurcation parameters
    a_cur = slider_a.val
    b_cur = slider_b.val

    #back parameters
    r_0 = 0.1
    m_0 = 0.2
    cVV_0 = slider_cVV.val
    e_0 = 0

    s_0 = slider_s.val
    n_0 = slider_n.val
    cSS_0 = 0.1

    u_0 = slider_u.val
    p_0 = slider_p.val
    cLL_0 = 0.1

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

E6_button.on_clicked(E6)

def E7(val):

    #bifurcation parameters
    a_cur = 0.25
    b_cur = 0.5

    #back parameters
    r_0 = slider_r.val
    m_0 = slider_m.val
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = slider_s.val
    n_0 = slider_n.val
    cSS_0 = slider_cSS.val

    u_0 = slider_u.val
    p_0 = slider_p.val
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

E7_button.on_clicked(E7)

def update2(val):

    # global slider_a.val, slider_b.val
    # #bifurcation parameters
    # slider_a.val = slider_anb.val
    # slider_b.val = slider_anb.val
    a_cur = slider_anb.val
    b_cur = slider_anb.val

    #back parameters
    r_0 = slider_r.val
    m_0 = slider_m.val
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = slider_s.val
    n_0 = slider_n.val
    cSS_0 = slider_cSS.val

    u_0 = slider_u.val
    p_0 = slider_p.val
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]


    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]

    x = odeint(dtChange, VSL_0, t, args=(a_cur, b_cur))

    V_plt.set_ydata(x[:,0])
    S_plt.set_ydata(x[:,1])
    L_plt.set_ydata(x[:,2])

    sim.canvas.draw()

def resetBackParam(val):
    slider_r.reset()
    slider_m.reset()
    slider_cVV.reset()
    slider_e.reset()

    slider_s.reset()
    slider_n.reset()
    slider_cSS.reset()

    slider_u.reset()
    slider_p.reset()
    slider_cLL.reset()

    update(None)

def resetBifurcation(val):
    slider_a.reset()
    slider_b.reset()

    update(None)

def resetPopulation(val):
    slider_V.reset()
    slider_S.reset()
    slider_L.reset()

    update(None)

#animated real time simulation figure

realTimeSim = plt.figure()
RTS = realTimeSim.subplots()
plt.subplots_adjust(bottom = 0.2)


t_sim, V_sim, S_sim, L_sim = [0], [VSL_0[0]], [VSL_0[1]], [VSL_0[2]]

def reset_simulation(val):
    global t_sim, V_sim, S_sim, L_sim, VSL_0
    VSL_0 = [slider_V.val, slider_S.val, slider_L.val]
    t_sim, V_sim, S_sim, L_sim = [0], [VSL_0[0]], [VSL_0[1]], [VSL_0[2]]


reset_RTsim = realTimeSim.add_axes([0.6, 0.05, 0.3, 0.04])
reset_simPop = Button(reset_RTsim, 'Reset Simulation', color='gold',hovercolor='skyblue')

def simulate(val):

    global t_sim, V_sim, S_sim, L_sim

    #bifurcation parameters
    a_cur = slider_a.val
    b_cur = slider_b.val

    #back parameters
    r_0 = slider_r.val
    m_0 = slider_m.val
    cVV_0 = slider_cVV.val
    e_0 = slider_e.val

    s_0 = slider_s.val
    n_0 = slider_n.val
    cSS_0 = slider_cSS.val

    u_0 = slider_u.val
    p_0 = slider_p.val
    cLL_0 = slider_cLL.val

    #redefine the function
    def dtChange(VSL_vector, t, a, b):
        #fox param
        r = r_0
        m = m_0
        cVV = cVV_0
        e = e_0
        #cottontail param
        s = s_0
        n = n_0
        cSS = cSS_0
        #hare param
        u = u_0
        p = p_0
        cLL = cLL_0

        V = VSL_vector[0]
        S = VSL_vector[1]
        L = VSL_vector[2]

        dVdt = V*(r - cVV*V - m + e*a*S + e*b*L)
        dSdt = S*(s - cSS*S - n - a*V)
        dLdt = L*(u - cLL*L - p - b*V)

        return [dVdt, dSdt, dLdt]

    VSL_sim0 = [V_sim[-1], S_sim[-1], L_sim[-1]]
    t_iter = np.linspace(t_sim[-1], t_sim[-1] + 3, 100)
    x = odeint(dtChange, VSL_sim0, t_iter, args=(a_cur, b_cur))

    t_sim = np.concatenate((t_sim[100:], t_iter), axis = None)
    V_sim = np.concatenate((V_sim[100:], x[:,0]), axis = None)
    S_sim = np.concatenate((V_sim[100:], x[:,1]), axis = None)
    L_sim = np.concatenate((V_sim[100:], x[:,2]), axis = None)
    
    RTS.cla()
    RTS.plot(t_sim, V_sim, 'r--')
    RTS.plot(t_sim, S_sim, 'blue')
    RTS.plot(t_sim, L_sim, 'black')
    RTS.legend(["Fox", "Cottontail", "Hare"])
    RTS.set_ylabel('Population')
    RTS.set_xlabel('time')

animation = FuncAnimation(realTimeSim, simulate, interval=200)

#actions

slider_a.on_changed(update)
slider_b.on_changed(update)
slider_anb.on_changed(update2)

slider_r.on_changed(update)
slider_m.on_changed(update)
slider_cVV.on_changed(update)
slider_e.on_changed(update)

slider_s.on_changed(update)
slider_n.on_changed(update)
slider_cSS.on_changed(update)

slider_u.on_changed(update)
slider_p.on_changed(update)
slider_cLL.on_changed(update)

slider_V.on_changed(update)
slider_S.on_changed(update)
slider_L.on_changed(update)

reset_back_param.on_clicked(resetBackParam)
reset_bifurcation.on_clicked(resetBifurcation)
reset_population.on_clicked(resetPopulation)
reset_simPop.on_clicked(reset_simulation)

    

plt.show()

