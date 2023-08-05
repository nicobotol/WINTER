function [z_sol1] = plane(A, B, C, xP, yP)
  %% Solve the system for finding the plane passing for 3 points in the space and then evaluate the z in the (vecotr of points) (xP, yP)
  
    syms z x y
  
    mat = [
      x-A(1) y-A(2) z-A(3); 
      B(1)-A(1) B(2)-A(2) B(3)-A(3); 
      C(1)-A(1) C(2)-A(2) C(3)-A(3)  
    ];
  
    det_mat = det(mat);
    zP = solve(det_mat == 0, z);
  
    z_sol = subs(zP, x = xP);
    z_sol1 = subs(z_sol, y = yP);
    z_sol1 = eval(z_sol1);
  end