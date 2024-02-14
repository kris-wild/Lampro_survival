## Steps for modelling

# Step 1: Get microclimate model for each lat and long
# Step 2: Ectotherm model for a mother based on the mircoclimate at that GPS point to simulate mother lifespane; age at maturity and reproductive events
# Step 3: Run DEB model. Identify the reproductive events
# Step 4: Extract how much energy / mass is put towards each reproductive event (Function for Step 3/4)
# Step 5: Extract the annual fecundity / number eggs produced
# Step 6: Now calculate the mass for a single egg based on Steps 4/5. This gives us the mass of a single / average egg
# Step 7: Using the daye of reproductive events, get soil / egg temperatures (10-15cm, on average, below ground) for 35 days AFTER each reproduction event.
# Step 8: Using temperature and egg mass from the specific site, use these to predict egg survival probaility based on experimental data
# Step 9: Repeat this across distribution

## Ultimate data
Lat, Long, ClimateScenerio, Egg mass, Egg temp/ Nest Temp, P surv (egg), sd(P surv), L95% (p-surv), U95% (p-surv)

# Step 10: Repeat this under climate change ERA projections
## Ultimate data
Lat, Long, ClimateScenerio, Egg mass, Egg temp/ Nest Temp, P surv (egg), sd(P surv), L95% (p-surv), U95% (p-surv)

##### In silico DEB experiment

1) Run ecotherm model, z constant and same as we have now, at two constant temps, 23 and 28C.
	- What is the body size at "hatch"/maturity?
	- What is the hatching time?

2) Run ecotherm model, z constant but Eo very low and high at two constant temps, 23 and 28C.
	- What is the MINIMUM Eo for development? Here you probably want to add in systematic variation to Eo (Eo varies from 100, 400, 600 you get idea.)
	- What is the body size at "hatch"/maturity?
	- What is the hatching time?

3) Run ecotherm model, z with rnorm(0,0.01) and same as we have now, at two constant temps, 23 and 28C.
	- What is the distribution of body size at "hatch"/maturity?
	- What is the distribution hatching time?

4) Run ecotherm model, z with rnorm(0,0.01) and same as we have now, at two constant temps, 23 and 28C with 15% less energy in Eo.
	- What is the distribution of body size at "hatch"/maturity?
	- What is the distribution hatching time?
