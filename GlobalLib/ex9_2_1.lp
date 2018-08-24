Minimize
  - 10 x(9) + 4 x(10) + [ 2 x(9) ^2 + 8 x(10) ^2 ] / 2 
Subject To
 con1: - 1.5 x(9) + 2 x(10) + x(5) - 0.5 x(6) + x(7) - x(8) = 2
 con2: - x(10) + x(4) = 0
 con3: x(9) + x(10) + x(3) = 7
 con4: x(9) - 0.5 x(10) + x(2) = 4
 con5: - 3 x(9) + x(10) + x(1) = -3
 QCP_row_for_con6: [ x(8) * x(4) ] = 0
 QCP_row_for_con7: [ x(7) * x(3) ] = 0
 QCP_row_for_con8: [ x(6) * x(2) ] = 0
 QCP_row_for_con9: [ x(5) * x(1) ] = 0
Bounds
 x(5) <= 20
 x(6) <= 20
 x(7) <= 20
 x(8) <= 20
 x(4) <= 20
 x(3) <= 20
 x(2) <= 20
 x(1) <= 20
End
