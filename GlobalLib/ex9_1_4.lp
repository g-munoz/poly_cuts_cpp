Minimize
  x(9) - 4 x(10)
Subject To
 con1: x(5) + 5 x(6) - 3 x(7) - x(8) = -1
 con2: - x(10) + x(4) = 0
 con3: 2 x(9) - 3 x(10) + x(3) = -4
 con4: 2 x(9) + 5 x(10) + x(2) = 108
 con5: - 2 x(9) + x(10) + x(1) = 0
 QCP_row_for_con6: [ x(5) * x(1) ] = 0
 QCP_row_for_con7: [ x(6) * x(2) ] = 0
 QCP_row_for_con8: [ x(8) * x(4) ] = 0
 QCP_row_for_con9: [ x(7) * x(3) ] = 0
Bounds
 x(5) <= 200
 x(6) <= 200
 x(7) <= 200
 x(8) <= 200
 x(4) <= 200
 x(3) <= 200
 x(2) <= 200
 x(1) <= 200
End
