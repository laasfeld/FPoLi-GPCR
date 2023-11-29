function scratchpad()
    RL = 3
    RRL = 3
    RR = 2
    RLRL = 1
    RRC = 1
    RCRC = 1
    RCRL = 1
    RL =1
    RC = 1
    R = 1
    L = 1
    nRR = 1
    
R1 = kf1*L*n*R-(((((kr1*RL)))))  ;                                                                                                                    
    R2 = kf2*C*n*R-(((((kr2*RC)))))   ;                                                                                                                    
    R3 = kf3*L*N-(((((kr3*NL)))))      ;                                                                                                                  
    R4 = kf4*L*NBV-(((((kr4*NBVL)))))   ;                                                                                                                  
    R5 = kf5*R*R*n*n-kr5*RR*nRR          ;                                                                                                                 
    R6 = kf6*RR*nRR*L-kr6*RRL             ;                                                                                                                
    R7 = kf7*RRL*L-kr7*RLRL               ;                                                                                                               
    R8 = kf8*RR*nRR*C-kr8*RRC               ;                                                                                                              
    R9 = kf9*RRC*C-kr9*RCRC                  ;                                                                                                           
    R10 = kf10*RRL*C-kr10*RCRL                ;                                                                                                            
    R11 = kf11*RRC*L-kr11*RCRL                 ;                                                                                                          
    R12 = kf12*RL*RL-kr12*RLRL                  ;                                                                                                          
    R13 = kf13*RC*RC-kr13*RCRC                  ;                                                                                                         
    R14 = kf14*RC*RL-kr14*RCRL ;
end