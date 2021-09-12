import scipy as sp
import pylab as plt
import math
from scipy.integrate import odeint

# Derived from HH original experimental data
def alpha_n(delta_voltage):
    V = float(delta_voltage)
    return (0.01 * (-10 -V))/(1 - math.exp((V + 10) / 10))

# Derived from HH original experimental data
def beta_n(delta_voltage):
    V = float(delta_voltage)
    return 0.125 * math.exp(V/80)

# Derived from HH original experimental data
def alpha_m(delta_voltage):
    V = float(delta_voltage)
    return (0.1 * (-25 - V))/(1 - math.exp((V + 25) / 10))

# Derived from HH original experimental data
def beta_m(delta_voltage):
    V = float(delta_voltage)
    return 4 * math.exp(V/18)

# Derived from HH original experimental data
def alpha_h(delta_voltage):
    V = float(delta_voltage)
    return 0.07 * math.exp(V/20)    

# Derived from HH original experimental data
def beta_h(delta_voltage):
    V = float(delta_voltage)
    return 1/(math.exp((V + 30) / 10) + 1)

class HodgkinHuxley():

    """Full Hodgkin-Huxley Model implemented in Python"""

    resting_potential = -65

    C_m  =   1.0
    """membrane capacitance, in uF/cm^2"""

    g_Na = 77.0
    """Sodium (Na) maximum conductances, in mS/cm^2"""

    g_K  =  24.0
    """Postassium (K) maximum conductances, in mS/cm^2"""

    g_L  =   0.3
    """Leak maximum conductances, in mS/cm^2"""

    E_Na =  50.0
    """Sodium (Na) Nernst reversal potentials, in mV"""

    E_K  = -77.0
    """Postassium (K) Nernst reversal potentials, in mV"""

    E_L  = -54.387
    """Leak Nernst reversal potentials, in mV"""


    resting_n = alpha_n(0) / (alpha_n(0) + beta_n(0)) 
    resting_m = alpha_m(0) / (alpha_m(0) + beta_m(0)) 
    resting_h = alpha_h(0) / (alpha_h(0) + beta_h(0)) 


    resting_g_k = g_K * resting_n ** 4
    resting_g_na = g_Na * (resting_m ** 3) * resting_h

    t = sp.arange(0.0, 250.0, 0.01)
    """ The time to integrate over """

    # def alpha_m(self, V):
    #     """Channel gating kinetics. Functions of membrane voltage"""
    #     return 0.1*(V+40.0)/(1.0 - sp.exp(-(V+40.0) / 10.0))

    # def beta_m(self, V):
    #     """Channel gating kinetics. Functions of membrane voltage"""
    #     return 4.0*sp.exp(-(V+65.0) / 18.0)

    # def alpha_h(self, V):
    #     """Channel gating kinetics. Functions of membrane voltage"""
    #     return 0.07*sp.exp(-(V+65.0) / 20.0)

    # def beta_h(self, V):
    #     """Channel gating kinetics. Functions of membrane voltage"""
    #     return 1.0/(1.0 + sp.exp(-(V+35.0) / 10.0))

    # def alpha_n(self, V):
    #     """Channel gating kinetics. Functions of membrane voltage"""
    #     return 0.01*(V+55.0)/(1.0 - sp.exp(-(V+55.0) / 10.0))

    # def beta_n(self, V):
    #     """Channel gating kinetics. Functions of membrane voltage"""
    #     return 0.125*sp.exp(-(V+65) / 80.0)


    def I_Na(self, V, m, h):
        """
        Membrane current (in uA/cm^2)
        Sodium (Na = element name)

        |  :param V:
        |  :param m:
        |  :param h:
        |  :return:
        """
        return self.g_Na * m**3 * h * (V - self.E_Na)

    def I_K(self, V, n):
        """
        Membrane current (in uA/cm^2)
        Potassium (K = element name)

        |  :param V:
        |  :param h:
        |  :return:
        """
        return self.g_K  * n**4 * (V - self.E_K)
    #  Leak
    def I_L(self, V):
        """
        Membrane current (in uA/cm^2)
        Leak

        |  :param V:
        |  :param h:
        |  :return:
        """
        return self.g_L * (V - self.E_L)

    def I_inj(self, t):
        """
        External Current

        |  :param t: time
        |  :return: step up to 10 uA/cm^2 at t>100
        |           step down to 0 uA/cm^2 at t>200
        |           step up to 35 uA/cm^2 at t>300
        |           step down to 0 uA/cm^2 at t>400
        """
        return 30*(t>48) - 30*(t>48.8) + 30*(t>148) - 30*(t>180) 


    #g = gt * n**4
    #n = (g/gt)^(1/4)
    # dgdt = 
    # dn\dt = d(g/gt^(1/4))\dt
    # dndt = (al)
    @staticmethod
    def dALLdt(X, t, self):
        """
        Integrate

        |  :param X:
        |  :param t:
        |  :return: calculate membrane potential & activation variables
        """
        V, m, h, n, gk, gna = X
        print("\n###################")
        print("Calculating dv/dt: current_v(%f):" % V)
        print("I_na current(%f), I_ka current(%f) I_l current(%f)" % (
                - self.I_Na(V, m, h),
                - self.I_K(V, n),
                - self.I_L(V)))

        gating_V = - (V - self.resting_potential)
        dVdt = (self.I_inj(t) - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_L(V)) / self.C_m
        dmdt = alpha_m(gating_V)*(1.0-m) - beta_m(gating_V)*m
        dhdt = alpha_h(gating_V)*(1.0-h) - beta_h(gating_V)*h
        dndt = alpha_n(gating_V)*(1.0-n) - beta_n(gating_V)*n

        print("dv/dt(%f): " %  dVdt)
        print("Calculating dn/dt: current_v(%f), current_n(%f):" % (V, n))
        print("dn/dt(%f):" % (dndt))
        print("Calculating dm/dt: current_v(%f), current_m(%f):" % (V, m))
        print("dm/dt(%f):" % (dmdt))
        print("Calculating dh/dt: current_v(%f), current_h(%f):" % (V, h))
        print("dh/dt(%f):" % (dhdt))

        dgkdt = (4 * (n ** 3) * dndt) * self.g_K
        dgnadt = self.g_Na * ((h * 3 * (m**2) * dmdt) + (m ** 3 * dhdt))

        # d[gna]/dt = gnamax(h * 3m ** 2 * dm/dt + m ** 3 * dh/dt)


        return dVdt, dmdt, dhdt, dndt, dgkdt, dgnadt 

    def Main(self):
        """
        Main demo for the Hodgkin Huxley neuron model
        """

        X = odeint(self.dALLdt, [-65, self.resting_m, self.resting_h, self.resting_n, self.resting_g_k, self.resting_g_na], self.t, args=(self,))
        V = X[:,0]
        m = X[:,1]
        h = X[:,2]
        n = X[:,3]
        gk = X[:,4]
        gNa = X[:,5]
        total = gk + gNa

        ina = self.I_Na(V, m, h)
        ik = self.I_K(V, n)
        il = self.I_L(V)

        plt.figure()

        plt.subplot(3,1,1)
        plt.title('Hodgkin-Huxley Neuron')
        plt.plot(self.t, V, 'k')
        plt.ylabel('V (mV)')
        plt.xlabel('t (ms)')
        #plt.set_xticks(0, 250, 20)
        plt.grid(color="white")
        ax = plt.gca()
        ax .set_facecolor("#ffd4d4")

        plt.subplot(3,1,2)
        plt.plot(self.t, gk, '#ffc800', label='$g_{K}$')
        plt.plot(self.t, gNa, 'r', label='$g_{Na}$')
        #plt.plot(self.t, total, 'black', label='$g_{T}$')
        plt.ylabel('Conductance ($S/cm^2$)')
        plt.xlabel('t (ms)')
        plt.legend()
        plt.grid(color="white")
        ax = plt.gca()
        ax .set_facecolor("#e8daeb")

        # plt.subplot(4,1,3)
        # plt.plot(self.t, ina, 'r', label='$I_{Na}$')
        # plt.plot(self.t, ik, 'b', label='$I_{K}$')
        # plt.plot(self.t, il, 'y', label='$I_{L}$')
        # plt.ylabel('Current')
        # plt.xlabel('t (ms)')
        # plt.legend()

        # plt.subplot(4,1,4)
        # plt.plot(self.t, m, 'r', label='m')
        # plt.plot(self.t, h, 'g', label='h')
        # plt.plot(self.t, n, 'b', label='n')
        # plt.ylabel('Gating Value (Probability)')
        # plt.legend()


        plt.subplot(3,1,3)
        i_inj_values = [self.I_inj(t) for t in self.t]
        plt.plot(self.t, i_inj_values, 'k')
        plt.xlabel('t (ms)')
        plt.ylabel('Injected current ($\\mu{A}/cm^2$)')
        plt.ylim(-1, 40)
        plt.grid(color="white")
        ax = plt.gca()
        ax .set_facecolor("#f0e2bd")
        plt.show()

if __name__ == '__main__':
    runner = HodgkinHuxley()
    runner.Main()