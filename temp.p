### SETUP DONE ###
pos_capacity = (c.POS_CSN_MAX - c.POS_CSN_INITIAL) * c.POS_ELEC_THICKNESS * (1-c.POS_ELEC_POROSITY) 
neg_capacity = (c.NEG_CSN_INITIAL - c.NEG_CSN_MIN) * c.NEG_ELEC_THICKNESS * (1-c.NEG_ELEC_POROSITY)

capacity = min(pos_capacity, neg_capacity) 
if (capacity == pos_capacity):
    print("Pos electrode parameters LIMITING")
else:
    print("Neg electrode parameters LIMITING")

capacity *= (c.F / 3600) # conversion into Ah
calc_current = (capacity / c.RUNTIME_HOURS)

seconds = c.RUNTIME_HOURS * 3600
time_steps = np.linspace(0, seconds, 250)
print(f"Evaluating @ {len(time_steps)} timesteps")
print(f"Discharging @ {calc_current:.3f} A/m2; Runtime: {seconds} seconds")