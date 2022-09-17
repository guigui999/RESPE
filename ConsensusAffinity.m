function [Z] = ConsensusAffinity(X,Z_ini,LM_ini,d,lambda1,lambda2,lambda3,lambda4,max_iter,rho,miu, max_miu)
         
        for i = 1:length(X)
                    
            % ---------- Initilization -------- %
            [m,n] = size(X{i});
            
            Y1{i} = zeros(n,n);
            Y2{i} = zeros(n,n);
            Y3{i} = zeros(d,d);           
            Y4{i} = zeros(m,n);
            Y5{i} = zeros(m,n);
            Y6{i} = zeros(d,n);

            
            R1{i} = ones(m,d);
            R2{i} = ones(m,d);
                     
            P1{i} = ones(m,d);
            P2{i} = ones(m,d);
            A{i} =  ones(d,d);
            
            E1{i}  = zeros(m,n);
            E2{i}  = zeros(m,n);
            E3{i}  = zeros(d,n);          
        end      

    for iter = 1:max_iter
           if iter == 1
                J = Z_ini;
                M = Z_ini;
                Z = Z_ini;
                LM = LM_ini;  
                AA = A;
                clear Z_ini  LM_ini
           end

            for i = 1:length(X)
                            
                [m,n] = size(X{i});
                
                % -------- Update P1 --------- %          
                B1 = X{i} * LM{i} * X{i}';
                C1 = X{i} - E1{i} +  Y4{i}/miu;
                G1 = 2 * lambda2 * B1 + 2* lambda3 * eye(m,m) + miu* X{i} * X{i}';   
                H1 = 2 * lambda3 * P2{i} * A{i} + miu * X{i}* C1'* R1{i};
                P1{i} = inv(G1) *  H1;
                clear B1 C1 G1 H1
                
                % -------- Update P2 --------- %
                B2 = X{i} - E2{i} +  Y5{i}/miu;
                C2 = X{i} * Z{i}; 
                D2 =  X{i} - X{i} * Z{i};
                F2 = E3{i} -  Y6{i}/miu;
                G2 = miu * (C2 * C2' + D2 * D2' )+ eye(m,m)* 1e-6;
                H2 = 2 * lambda3 * (A{i}* A{i}') ;
                J2 = 2 * lambda3 * P1{i}* A{i}' + miu * C2*  B2' *  R2{i}+ miu * D2 * F2';  
                
                P2{i} = lyap(G2, H2, -J2);
                clear B2 C2 D2 F2 G2 H2  J2
                
                % -------- Update A --------- %  
                temp_0 = 2 * lambda3 * P2{i}'*P2{i} + miu * eye(d,d);
                temp_1 = 2 * lambda3 * P2{i}'*P1{i} + miu * AA{i} + Y3{i};
                A{i}= inv(temp_0) * temp_1;   
                clear  temp_0 temp_1
                
                 % -------- Update AA --------- %  
                eps1 = lambda4/miu ;
                temp_00 = A{i} - Y3{i}/miu;
                AA{i} = max(0,temp_00 - eps1) + min(0, temp_00 + eps1);
                clear  eps1 temp_00
                
                % -------- Update R1 --------- %                          
                R1_0 = X{i} - E1{i} +  Y4{i}/miu;
                G =  R1_0  * X{i}' * P1{i};       
                [U1,~,V1] = svd(G,'econ');
                R1{i} = U1 * V1';
                clear  R1_0 G  U1 V1

              % -------- Update R2 --------- %                             
                R2_0 = X{i} - E2{i} +  Y5{i}/miu;
                G   =  R2_0 *Z{i}'* X{i}' * P2{i};  
                [U2,~,V2] = svd(G,'econ');
                R2{i} = U2 * V2';
                clear R2_0  G U2 V2     
                
              % -------- Update J --------- %                             
                J{i} = softth((Z{i}-Y1{i}/miu) + eye(n) * 1e-8,1/miu );% J = softth((Z-Y4/mu)+eye(N)*1e-8,lambda/mu);
             
               % -------- Update M --------- %          
                distX{i} = L2_distance_1(P1{i}'*X{i},P1{i}'*X{i});
                dist  = 0.5 * lambda2 * distX{i};
                M{i}  = Z{i} - (Y2{i}+dist)/miu;
                for ic = 1:n
                    idx    = 1:n;
                    idx(ic) = [];
                    M{i}(ic,idx) = EProjSimplex_new(M{i}(ic,idx));          % 
                end  
                
                M00 = (M{i}+M{i}')/2;
                D00 = diag(sum(M00));
                LM{i} = D00 - M00;             
                clear M00 D00 distX dist ic idx
               
                % -------- Update Z --------- %
                 B4 = X{i} - E2{i} +  Y5{i}/miu;
                 A4 = P2{i}' *  X{i}- E3{i}+ Y6{i}/miu;
                 C4 = 2 * eye (n) + 2 *  X{i}' * P2{i}* P2{i}' *  X{i}; 
                 D4 = J{i} + M{i};
                 F4 = (Y1{i} + Y2{i})/miu;
                 G4 =  X{i}' * P2{i}* R2{i}' * B4;
                 H4 = X{i}' * P2{i}* A4;
                 J4 = D4 + F4 + G4+ H4;
                 Z{i} = inv(C4)* J4;   %Z{i} = inv(G4)* L4;
                 clear A4 B4 C4 G4 D4 F4 G4 H4 J4
                
                % ------- Update E ---------- %
                 B5 = X{i} - R1{i}*P1{i}' * X{i} + Y4{i}/miu;
                 C5 = X{i} - R2{i}*P2{i}' * X{i}* Z{i} + Y5{i}/miu;
                 D5 =  P2{i}'*X{i} - P2{i}' * X{i}* Z{i} + Y6{i}/miu;
                 G5 = [B5; C5; D5];
                 [E] = solve_l1l2(G5,lambda1/miu);
                 E1{i} = E(1:m,:); 
                 E2{i} = E(1+m:2*m,:);
                 E3{i} = E(1+2*m:2*m+d,:);
                 clear B5 C5 D5 G5 E
                                                        
                % -------- Update Y1 Y2 Y3 Y4 Y5 Y6-------- %           
                Y1{i} = Y1{i} + miu*(J{i}-Z{i});
                Y2{i} = Y2{i} + miu*(M{i}-Z{i});
                Y3{i} = Y3{i} + miu*(AA{i}- A{i});      
                Y4{i} = Y4{i} + miu*( X{i} - R1{i} * P1{i}' * X{i} - E1{i});
                Y5{i} = Y5{i} + miu*( X{i} - R2{i} * P2{i}' * X{i}* Z{i} - E2{i});
                Y6{i} = Y6{i} + miu*( P2{i}'*X{i} - P2{i}' * X{i}* Z{i} - E3{i});
            end
                                      
               LL1 = sparse(1,length(X));
               LL2 = LL1;
               LL3 = LL1;
               LL4 = LL1;
               LL5 = LL1;
               LL6 = LL1;
                            
            for i = 1:length(X)
                LL1(i) = norm(J{i}-Z{i},inf);
                LL2(i) = norm(M{i}-Z{i},inf);
                LL3(i) = norm(AA{i}-A{i},inf);
                LL4(i) = norm( X{i} - R1{i}*P1{i}' * X{i} - E1{i},inf);
                LL5(i) = norm( X{i} - R2{i}*P2{i}' * X{i}* Z{i} - E2{i},inf);   
                LL6(i) = norm(P2{i}'*X{i} - P2{i}' * X{i}* Z{i} - E3{i},inf);  
            end
                           
            % ---------- miu && convergence ------------- %
            miu = min(rho*miu,max_miu);
           if ( max(LL1) < 1e-6 && max(LL2) < 1e-6 && max(LL3) < 1e-6 && max(LL4) < 1e-6 && max(LL5)< 1e-6 && max(LL6)< 1e-6)
           iter
               break;
           end    
              
    end
   clear  Y1 Y2 Y3 Y4 Y5 R1 R2 P1 P2 A AA E1 E2 E3 J M LM  LL1 LL2 LL3 LL4 LL5 distX dist d miu m n miu 
end