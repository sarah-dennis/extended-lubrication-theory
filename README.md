# EXTENDED LUBRICATION THEORY
Code for, 'Comparison of Extended Lubrication Theories for Stokes Flow' (Sarah Dennis & Thomas Fai, 2026)

This project includes finite difference solvers for the Reynolds equation, the Stokes equation,
and various models in extended/perturbed lubrication theory.

Details on the extended lubrication theory models (VA-ELT and T.G.-ELT), and the 
perturbed lubrication theory models (epsilon^2-PLT and epsilon^4-PLT) can be found in the paper.

#------------------------------------------------------------------------------
# REYNOLDS, ELT & PLT SOLVERS
#------------------------------------------------------------------------------

To use the Reynolds, PLT, or ELT solvers, use file reyn_run.py

1. Choose an Example and set the geometry parameters (args)
    
    A.) backward facing step

          #--------------------------------------------------------------------
            Example = examples.BFS
            h_inlet = 2
            h_outlet = 1
            l_inlet = 8
            L_total = 16
            args =  [h_inlet, h_outlet, l_inlet, L_total]
          #--------------------------------------------------------------------
    
    B.) triangular textured slider
    
          #--------------------------------------------------------------------
            Example = examples.TriSlider
            h_in = 1
            h_vertex = 2
            h_out = h_in
            l_in = 7
            l_out = 7
            l_a = 1.25
            l_b = 0.75
            args =  [h_in, h_vertex, h_out, l_in, l_a, l_b, l_out]
          #--------------------------------------------------------------------
    
    C.) logistic step
    
          #--------------------------------------------------------------------
            Example = examples.Logistic
            delta = 8     # max slope: delta*(h_in-h_out)/4
            h_in = 1       
            h_out = 2    
            L_total = 16     
            args = [ h_in, h_out, L_total, delta]
          #--------------------------------------------------------------------

    --> See reyn_examples.py, reyn_heights.py, domain.py
    

2. Set the boundary conditions 

    - U: lower surface velocity
    - Q: flux 
    
          #--------------------------------------------------------------------    
            U = 0
            Q = 1
            BC = bc.Mixed(U,Q_)
          #--------------------------------------------------------------------   
    
    --> See boundary.py

3. Initialize the Solver and grid scale (N)

    - Solver: class with Reynolds, ELT, PLT solve methods
    - N: grid scale (dx = 1/N)
    
          #--------------------------------------------------------------------    
            Solver = control.Reynolds_Solver(Example, BC, args)
            N = 160
          #--------------------------------------------------------------------    
          
    --> see reyn_control.py

4. Run the finite difference method(s)

    - Each solve method will plot the pressure and velocity solution
    - The bool-flags plots_on, zoom_on are set at the top of reyn_run.py
   
    A.) Reynolds solver
    
          #--------------------------------------------------------------------    
            Solver.fd_solve(N, plot=plots_on, zoom=zoom_on)
          #--------------------------------------------------------------------
    
      --> see reyn_finDiff.py, reyn_velocity.py   
           
    B.) T.G.-ELT solver 
    
          #--------------------------------------------------------------------   
            Solver.fd_TG_ELT_solve(N, plot=plots_on, zoom=zoom_on)
          #--------------------------------------------------------------------   
        
      --> see reyn_pressure_ELT.py, reyn_velocity.py
        
    C.) VA-ELT solver
    
          #--------------------------------------------------------------------   
            Solver.fd_VA_ELT_solve(N, plot=plots_on, zoom=zoom_on)
          #--------------------------------------------------------------------  
          
      --> see reyn_pressure_ELT.py, reyn_velocity_ELT.py
          
    D.) PLT solver
    
          #--------------------------------------------------------------------   
            Solver.fd_pert_solve(N, order=2,  plot=plots_on,  zoom=zoom_on)
            
            Solver.fd_pert_solve(N, order=4,  plot=plots_on, zoom=zoom_on)
          #--------------------------------------------------------------------   

      --> see reyn_perturbed.py
       
   --> see reyn_control.py,  graphics.py
 
#------------------------------------------------------------------------------
# STOKES SOLVER
#------------------------------------------------------------------------------

To use the Stokes solver, use file stokes_run.py

1. Choose an Example and set the geometry parameters (args)
    
    - Examples are set identically to reyn_run.py above
    
    A.) backward facing step (examples.BFS)
    
    B.) triangular textured slider
    
    C.) logistic step (examples.logistic)
    
    --> See stokes_examples.py, stokes_heights.py, domain.py
    

2. Set the boundary conditions 

    - Boundary conditions are set identically to reyn_run.py above

    - U: lower surface velocity
    - Q: flux 

    --> See boundary.py

3. Initialize the Solver and grid scale (N)

    - Solver: class with solve methods
    - N: grid scale (dx = 1/N)
    - Re: Reynolds number for stream-velocity Navier-Stokes equations (use Re=0 for Stokes solver)
    
          #--------------------------------------------------------------------    
            Solver = control.Stokes_Solver(Example, args, BC, Re=0) 
            
            N = 160
          #--------------------------------------------------------------------    
          
    --> see stokes_control.py

4. Run the iterative Stokes solver

    - The iterative Stokes solver computes the stream function and velocity, 
      in general the pressure is not computed
    - The iterative solver will run until the error in the stream function between
      iterations falls below 10^-8, or until the maximum number of iterations has been reached
    - The solver checks this error and saves the current solution every 500 iterations
    - The solver will also print the current error and iteration number to the console
    - The parameters {err_tol, max_iters, error_mod, write_mod} can be updated in stokes_control.py
    - solutions are written to ./stokes_examples/example_name/example_name_N

    A.) run a new solution at grid size N
    
          #--------------------------------------------------------------------    
            Solver.new_run(N)
          #--------------------------------------------------------------------    
    
    B.) load and run an existing solution at grid size N
    
          #--------------------------------------------------------------------    
            Solver.load_run(N)
          #--------------------------------------------------------------------      
          
    C.) run a new solution starting at grid size N,
        then interpolate the solution to grid size 2*N and run,
        repeat, ending with grid size (2^k)*N
    
          #--------------------------------------------------------------------
            k = 4    
            Solver.new_run_many(N, 2, k)
          #--------------------------------------------------------------------    
    
    D.) load and run an existing solution starting at grid size N,
        then interpolate the solution to grid size 2*N and run,
        repeat, ending with grid size (2^k)*N   

          #--------------------------------------------------------------------
            k = 4    
            Solver.load_run_new_many(N, 2, k)
          #--------------------------------------------------------------------    
        
    E.) load and run an existing solution starting at grid size N,
        then load and run an existing solution at grid size 2*N,
        repeat, ending with grid size (2^k)*N   

          #--------------------------------------------------------------------
            k = 4    
            Solver.load_run_many(N, 2, k)
          #--------------------------------------------------------------------    
    
    F.) load an existing solution at grid size N
        
        - this is the only method that computes the pressure and returns the solution

          #--------------------------------------------------------------------   
            p, u, v, stream = Solver.load(N)
          #--------------------------------------------------------------------    

    --> see stokes_control.py, stokes_solver.py, stokes_pressure.py, read_write.py

5. Plot the solution
    
    - the solution at grid size N must already exist (see step 4.)
    - the pressure is computed from the saved stream-velocity solution
    - the bool-flag zoom-on is set at the top of stokes_run.py
    
          #--------------------------------------------------------------------
            Solver.load_plot(N, zoom=zoom_on)
          #--------------------------------------------------------------------    
       
    --> see stokes_control.py,  graphics.py

 
#------------------------------------------------------------------------------
# SOLUTION COMPARISION
#------------------------------------------------------------------------------

To compare the velocity and pressure solutions from the Reynolds, ELT, PLT and 
Stokes methods, use the file compare.py

The 5 methods are run for a specific example with one geometric parameter varied

The l2 norm relative error in pressure and velocity is plotted

1. Choose an Example and set the geometric parameters

    - Example is set similarly to reyn_run.py and stokes_run.py above
    - Set test_args with the geometric parameter being varied 
    
    A.) Logistic step with varying delta
    
          #--------------------------------------------------------------------
            Reyn_Example = reyn_examples.Logistic
            Stokes_Example = stokes_examples.Logistic
            h_in = 2   
            h_out= 1   
            l_total = 16   #  length
            
            tests = [2, 3, 4, 6, 8, 16, 32]
            test_args = [[h_in , h_out, l_total, delta] for delta in tests]
            
            exstr = 'Logistic Step'
            label = '$\lambda$'
          #--------------------------------------------------------------------

    B.)  Triangular slider with varying H_vertex
    
          #--------------------------------------------------------------------
            Stokes_Example = stokes_examples.TriSlider
            
            h_in=1  
            h_out = 1 
            l_in = 7  
            l_out = 7  
            l_a = 1.25 
            l_b = 0.75  
            
            tests = [1/16, 1/8, 1/4, 1/2, 3/4, 5/4, 3/2, 7/4, 2]
            test_args = [[h_in, h_vertex, h_out, l_in, l_a, l_b, l_out] for h_vertex in tests]
            
            exstr = 'Triangular Slider'
            label = '$H_v$'
          #--------------------------------------------------------------------

2. Set the boundary conditions (BC) and grid scale (N)
    
    - Identically as in reyn_run.py and stokes_run.py

3. Ensure stokes solutions are precomputed

3. Run the script compare.py

    - plots the l2 norm relative error for pressure and for velocity 
    
    - adjust plot_2D in graphics.py to set y_lim, x_lim, and legend location as needed


