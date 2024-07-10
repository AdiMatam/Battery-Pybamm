import pybamm
import consts as cc

pybamm.set_logging_level("DEBUG")

#OCP_INIT = 4.027031539154782
NEG_OCP_INIT = 0.08352811644995728
POS_OCP_INIT = 4.08138601219583

class SingleParticle:
    def __init__(self, name: str, charge: int, iapp: pybamm.Variable, r_val: float):
        self.name = name
        self.domain = name + " dDomain"
        self.charge = charge

        self.c0 = pybamm.Parameter(name + " pInitial Concentration")
        self.L = pybamm.Parameter(name + " pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + " pElectrode Porosity")
        self.cmax = pybamm.Parameter(name + " pMax Concentration")

        self.iflag = pybamm.Parameter(name + " pCharge?")

        self.iapp = iapp

        self.i_s = pybamm.Variable(name + " vSide Current")
        self.i_int = pybamm.Variable(name + " vIntercalation Current")
        self.seiL = pybamm.Variable(name + " vSEI Thickness")
        self.sei0 = pybamm.Parameter(name + " pInitial SEI Length")

        self.phi = pybamm.Variable(name + " vPhi")
        self.csn = pybamm.Variable(name + " vConcentration", domain=self.domain)
        self.surf_csn = pybamm.surf(self.csn)

        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")
        self.r_val = r_val

    def u_func(self, sto):
        return pybamm.FunctionParameter(
            self.name + " fParticle Voltage", 
            {
                "Stoichiometry": sto
            }
        )

    def process_model(self, model: pybamm.BaseModel, clear=False):
        flux = cc.D * -pybamm.grad(self.csn)
        # dc/dt = d^2c/dr^2
        dcdt = -pybamm.div(flux)

        # solve the diffusion equation (del * del(c))
        model.rhs.update({
            self.csn: dcdt,
        })

        # self.charge = +1 for positive electrode
        #               -1 for negative electrode
        self.a_term = (3 * (1 - self.eps_n)) / self.r_val
        self.j = (self.charge * self.iapp) / (self.L * self.a_term)

        self.ocp = self.u_func(self.surf_csn / self.cmax)

        cathode_kint = 1.04e-11
        j0 = cc.F * cathode_kint * self.surf_csn**0.5 * (self.cmax - self.surf_csn)**0.5 
        x = cc.F / (2 * cc.R_GAS * cc.T) * (self.phi - self.ocp)

        model.algebraic.update({
            self.phi: j0 * 2 * pybamm.sinh(x) - self.j,
        })

        model.initial_conditions.update({
            self.csn: pybamm.x_average(self.c0),
            self.phi: POS_OCP_INIT
        }) 

        model.boundary_conditions.update({
            self.csn: {
                "left":  (0, "Neumann"),
                "right": (-self.j / (cc.F * cc.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })

        model.variables.update({
            self.csn.name: self.csn,
            self.phi.name: self.phi
        })


    def process_anodeNOSEI(self, model: pybamm.BaseModel, clear=False):
        flux = cc.D * -pybamm.grad(self.csn)
        # dc/dt = d^2c/dr^2
        dcdt = -pybamm.div(flux)

        # solve the diffusion equation (del * del(c))
        model.rhs.update({
            self.csn: dcdt,
        })

        # self.charge = +1 for positive electrode
        #               -1 for negative electrode
        self.a_term = (3 * (1 - self.eps_n)) / self.r_val
        self.j = (self.charge * self.iapp) / (self.L * self.a_term)

        self.ocp = self.u_func(self.surf_csn / self.cmax)

        anode_kint = 2.07e-11
        j0 = cc.F * anode_kint * self.surf_csn**0.5 * (self.cmax - self.surf_csn)**0.5 
        x = cc.F / (2 * cc.R_GAS * cc.T) * (self.phi - self.ocp)

        model.algebraic.update({
            self.phi: j0 * 2 * pybamm.sinh(x) - self.j,
        })

        model.initial_conditions.update({
            self.csn: pybamm.x_average(self.c0),
            self.phi: NEG_OCP_INIT
        }) 

        model.boundary_conditions.update({
            self.csn: {
                "left":  (0, "Neumann"),
                "right": (-self.j / (cc.F * cc.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })

        model.variables.update({
            self.csn.name: self.csn,
            self.phi.name: self.phi
        })

    def process_anode(self, model: pybamm.BaseModel, clear=False):
        flux = cc.D * -pybamm.grad(self.csn)
        # dc/dt = d^2c/dr^2
        dcdt = -pybamm.div(flux)

        # self.charge = +1 for positive electrode
        #               -1 for negative electrode
        self.a_term = (3 * (1 - self.eps_n)) / self.r_val
        self.j = (self.charge * self.iapp) / (self.L * self.a_term)

        KSEI = 5.0e-6
        M_SEI = 0.162
        RHO_SEI = 1690

        ## -- SEI START -- 
        dLdt = (-self.i_s / 2*cc.F) * (M_SEI / RHO_SEI)

        # solve the diffusion equation (del * del(c))
        model.rhs.update({
            self.csn: dcdt,
            self.seiL: dLdt
        })

        ## anode (SEI case)
        self.ocp = self.u_func(self.surf_csn / self.cmax)
        x = cc.F / (2 * cc.R_GAS * cc.T) * (self.phi - self.ocp - (self.seiL/KSEI)*self.j)

        kfs = 1.36e-12
        ## TODO: CHECK WITH PROF: e_sei * ce_sln (taken from paper)
        cec_init = 0.05 * 4.541
        is_rhs = self.iflag * -cc.F*kfs*cec_init * pybamm.exp( (-0.5*cc.F)/(cc.R_GAS*cc.T) * (self.phi - (self.seiL/KSEI)*self.j) ) 

        anode_kint = 2.07e-11
        j0 = cc.F * anode_kint * self.surf_csn**0.5 * (self.cmax - self.surf_csn)**0.5 
        model.algebraic.update({
            self.phi: j0 * 2 * pybamm.sinh(x) - self.i_int,
            self.i_s: is_rhs - self.i_s,
            self.i_int: self.i_int + self.i_s - self.j
        })

        model.initial_conditions.update({
            self.csn: pybamm.x_average(self.c0),
            self.phi: NEG_OCP_INIT,
            self.i_s: 0.5 * self.j,
            self.i_int: 0.5 * self.j,
            self.seiL: self.sei0,
        }) 

        model.boundary_conditions.update({
            self.csn: {
                "left":  (0, "Neumann"),
                "right": (-self.i_int / (cc.F * cc.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })

        model.variables.update({
            self.csn.name: self.csn,
            self.phi.name: self.phi,
            self.i_int.name: self.i_int,
            self.i_s.name: self.i_s,
            self.seiL.name: self.seiL,
            "Jterm": self.j
        })
    
    def process_geometry(self, geo: dict, clear=False):
        if clear:
            geo.clear()

        ## TODO: need another 'domain' for SEI
        geo.update({
            self.domain: {self.r: {"min": 0, "max": self.r_val}}
        })
    
    def process_parameters(self, all_params: dict, particle_params: dict, clear=False):
        if clear:
            all_params.clear()

        all_params.update(
            {key.name : value for key, value in particle_params.items()}
        )









if __name__ == '__main__':
    import params as p
    import numpy as np
    from ocp import NEG_OCP, POS_OCP

    I_TOTAL = +1.2027659291666666
    DISCRETE_PTS = 30
    TIME_PTS = 250

    geo = {}
    model = pybamm.BaseModel()

    i_total = pybamm.Parameter("Input Current") 
    parameters = {}

    a = SingleParticle("Anode", -1, i_total, p.PARTICLE_RADIUS.get_value())
    a.process_anodeNOSEI(model)
    a.process_geometry(geo)

    a.process_parameters(parameters, {
        i_total:            "[input]",
        a.c0:               "[input]",
        a.L:                p.NEG_ELEC_THICKNESS.get_value(),
        a.eps_n:            p.NEG_ELEC_POROSITY.get_value(),
        a.cmax:             p.NEG_CSN_MAX.get_value(),

        a.ocp:              p.NEG_OCP,
        #a.sei0:             "[input]",
        #a.iflag:            "[input]",
    })

    particles = [a]
    mesh = pybamm.Mesh(geo, 
        { d.domain: pybamm.Uniform1DSubMesh for d in particles },
        { d.r: DISCRETE_PTS for d in particles }
    )

    disc = pybamm.Discretisation(mesh, 
        { d.domain: pybamm.FiniteVolume() for d in particles }
    )

    param_ob = pybamm.ParameterValues(parameters)
    param_ob.process_model(model)
    param_ob.process_geometry(geo)

    disc.process_model(model)

    cycles = 1
    solver = pybamm.CasadiSolver(mode='safe', atol=1e-6, rtol=1e-5, extra_options_setup={"max_num_steps": 100000})
    HOURS = 10
    time_steps = np.linspace(0, 3600 * HOURS, TIME_PTS)
    total_time_steps = np.linspace(0, 3600 * HOURS * cycles, TIME_PTS * cycles)
    
    sign = -1
    last_conc = p.NEG_CSN_INITIAL.get_value()
    last_SEI = 5.e-9

    conc_array = np.empty((TIME_PTS*cycles,))
    volt_array = np.empty((TIME_PTS*cycles,))
    #int_array = np.empty((TIME_PTS*cycles,))
    #sei_array = np.empty((TIME_PTS*cycles,))
    #seil_array = np.empty((TIME_PTS*cycles,))

    for i in range(cycles):
        inps = {
            i_total.name: sign * I_TOTAL,
            a.c0.name: last_conc,
            #a.sei0.name: last_SEI,
            #a.iflag.name: 0 if (sign == -1) else 1
        }
        solution = solver.solve(model, time_steps, inputs=inps)

        arr = solution[a.csn.name].entries[-1]
        conc_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = arr
        volt_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = solution[a.phi.name].entries
        #int_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = solution[a.i_int.name].entries
        #sei_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = solution[a.i_s.name].entries
        #seil_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = solution[a.seiL.name].entries

        last_conc = arr[-1]
        #last_SEI = solution[a.seiL.name].entries[-1]
        sign *= -1
        print(f"Passed {i}")

    print(len(conc_array))
    print(len(total_time_steps))

    import pandas as pd
    df = pd.DataFrame({
        'Time': total_time_steps,
        'Concentration': conc_array,
        "Voltage": volt_array,
        #'I_INT': int_array,
        #'I_SEI': sei_array,
        #"Length": seil_array,
    })

    df.to_csv("ANODE_NOSEI_1.csv", index=False)

    quit()


    from matplotlib import pyplot as plt

    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Time (s)')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:red'
    ax1.set_ylabel('Concentration (mol / m2)', color=color)
    ax1.plot(total_time_steps, conc_array, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2.set_ylabel('Voltage (V)', color='tab:blue')  # we already handled the x-label with ax1
    ax2.plot(total_time_steps, volt_array, color='blue')
    ax2.tick_params(axis='y', labelcolor='tab:blue')

    plt.grid(linewidth=0.2)
    plt.tight_layout()
    plt.show()



        ## Butler volmer calculation
        # self.ocp = self.u_func(self.surf_csn / self.cmax)
        # self.overp = (2*c.RTF*pybamm.arcsinh(self.j / (2 * self.j0)))
        # self.phi = self.ocp + self.overp