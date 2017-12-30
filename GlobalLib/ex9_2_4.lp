Minimize
  - 2 x(5) - 2 x(6) + [ x(5) ^ 2 + x(6) ^ 2 ] / 2 
Subject To
 con1: - x(4) + x(7) = -1
 con2: x(5) + x(7) - x(3) = 0
 con3: - x(6) + x(2) = 0
 con4: - x(5) + x(1) = 0
 con5: x(5) + x(6) - x(8) = 0
 QCP_row_for_con6: [ x(4) * x(2) ] = 0
 QCP_row_for_con7: [ x(3) * x(1) ] = 0
Bounds
 x(4) <= 200
 x(7) free
 x(3) <= 200
 x(2) <= 200
 x(1) <= 200
End
