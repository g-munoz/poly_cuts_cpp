Minimize
  - 3 x(13) - 3 x(14) + 2 x(15) + 2 x(16)
Subject To
 con1: 2 x(14) - 2 x(16) + 2 x(8) - x(11) + x(12) = -40
 con2: 2 x(13) - 2 x(15) + 2 x(7) - x(9) + x(10) = -40
 con3: x(14) + x(6) = 20
 con4: - x(14) + x(5) = 10
 con5: x(13) + x(4) = 20
 con6: - x(13) + x(3) = 10
 con7: 2 x(14) - x(16) + x(2) = -10
 con8: 2 x(13) - x(15) + x(1) = -10
 con9: x(13) - 2 x(14) + x(15) + x(16) <= 40
 QCP_row_for_con10: [ x(7) * x(1) ] = 0
 QCP_row_for_con11: [ x(8) * x(2) ] = 0
 QCP_row_for_con12: [ x(9) * x(3) ] = 0
 QCP_row_for_con13: [ x(10) * x(4) ] = 0
 QCP_row_for_con14: [ x(12) * x(6) ] = 0
 QCP_row_for_con15: [ x(11) * x(5) ] = 0
Bounds
 x(13) free
 x(14) free
 x(15) <= 50
 x(16) <= 50
 x(8) <= 200
 x(11) <= 200
 x(12) <= 200
 x(7) <= 200
 x(9) <= 200
 x(10) <= 200
 x(6) <= 200
 x(5) <= 200
 x(4) <= 200
 x(3) <= 200
 x(2) <= 200
 x(1) <= 200
End
