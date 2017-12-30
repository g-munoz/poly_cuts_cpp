Minimize
  2 x(5) + 3 x(6) + [ - 8 x(5) * x(6) ] / 2 
Subject To
 con1: 4 x(5) - x(3) + x(4) = 1
 con2: x(6) + x(2) = 1
 con3: - x(6) + x(1) = 0
 QCP_row_for_con4: [ x(4) * x(2) ] = 0
 QCP_row_for_con5: [ x(3) * x(1) ] = 0
Bounds
 x(5) <= 1
 x(3) = 0
 x(4) = 0
 x(2) <= 20
 x(1) <= 20
End
