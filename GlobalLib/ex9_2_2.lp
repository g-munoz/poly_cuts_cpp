Minimize
  0 x(4) + 0 x(8) - 20 x(10) + [ 2 x(10) ^ 2 + 2 x(9) ^ 2 ] / 2 
Subject To
 con1: 4 x(10) + 2 x(9) + x(5) - x(6) + x(7) = 60
 con2: x(10) + x(3) = 20
 con3: - x(10) + x(2) = 0
 con4: x(10) + x(9) + x(1) = 20
 con5: x(10) - x(9) <= 0
 QCP_row_for_con6: [ x(4) * x(8) ] = 0
 QCP_row_for_con7: [ x(5) * x(1) ] = 0
 QCP_row_for_con8: [ x(7) * x(3) ] = 0
 QCP_row_for_con9: [ x(6) * x(2) ] = 0
Bounds
 x(4) <= 20
 x(8) <= 20
 x(9) <= 15
 x(5) <= 20
 x(6) <= 20
 x(7) <= 20
 x(3) <= 20
 x(2) <= 20
 x(1) <= 20
End
