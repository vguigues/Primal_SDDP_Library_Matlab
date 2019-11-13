

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T: number of stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size_x(t): number of x variables for stage t, t=1,..,T.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size_u(t): number of u variables for stage t, t=1,..,T.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size_b(t): number of components in b(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upper_x: upper_x{1,t} is the upper bound on x(t) for stage t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lower_x: lower_x{1,t} is the lower bound on x(t) for stage t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upper_u: upper_u{1,t} is the upper bound on u(t) for stage t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lower_u: lower_u{1,t} is the lower bound on u(t) for stage t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%costc{1,t}: cost vector c for stage t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%costd{1,t}: cost vector d for stage t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bs{1,t}{1,j}: rhs b_t for stage t scenario j.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fs{1,t}{1,j}: rhs f_t for stage t scenario j.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The technology matrix A_{t} for stage t is given by
%a_subi{1,t}, a_subj{1,t}, a_valij{1,t}:
%A_{t}(a_subi{1,t}(k),a_subj{1,t}(k))=a_valij{1,t}(k).
%Similarly for matrices B_{t}, C_{t} for stage t>=2. For instance B_t is
%given by b_subi{1,t-1}, b_subj{1,t-1}, b_valij{1,t-1}:
%B_{t}(b_subi{1,t-1}(k),b_subj{1,t-1}(k))= b_valij{1,t-1}(k).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%probabilities{1,t-1}(j) is the probability of scenario j for stage t=2,..,T.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nb_scenarios_rhs(t): number of scenarios for stage t=1,..,T.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nb_iter_max: maximal number of iterations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upper_bounds(k) is the upper bound for iteration k.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lower_bounds(k) is the lower bound for iteration k.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time(k) is the CPU time needed to solve iteration k of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lower_bounds,upper_bounds,time]=primal_sddp_control_rhs_eq(T,size_x,size_u,size_b,upper_x,lower_x,upper_u,lower_u,costc,costd,bs,probabilities,nb_scenarios_rhs,nb_iter_max,a_subi,a_subj,a_valij,b_subi,b_subj,b_valij,c_subi,c_subj,c_valij,big_m,talpha,tol)

lower_bounds=[];
upper_bounds=[];
time=[];

dynamic_subi=cell(1,T);
dynamic_subj=cell(1,T);
dynamic_valij=cell(1,T);

dynamic_buc=cell(1,T);
dynamic_blc=cell(1,T);

dynamic_blx=cell(1,T);
dynamic_bux=cell(1,T);

dynamic_c=cell(1,T);

%Initializations

dynamic_subi{1,1}=[a_subi{1,1},c_subi{1,1},size_b(1)+1];
dynamic_subj{1,1}=[a_subj{1,1},c_subj{1,1}+size_x(1),size_x(1)+size_u(1)+1];
dynamic_valij{1,1}=[a_valij{1,1},c_valij{1,1},1];
dynamic_blx{1,1}=[lower_x{1,1};lower_u{1,1};-inf];
dynamic_bux{1,1}=[upper_x{1,1};upper_u{1,1};inf];
dynamic_blc{1,1}=[bs{1,1}{1,1};-big_m];
dynamic_buc{1,1}=[bs{1,1}{1,1};inf];
dynamic_c{1,1}=[costc{1,1};costd{1,1};1];

for t=2:T-1
    dynamic_subi{1,t}=[a_subi{1,t},c_subi{1,t},size_b(t)+1];
    dynamic_subj{1,t}=[a_subj{1,t},c_subj{1,t}+size_x(t),size_x(t)+size_u(t)+1];
    dynamic_valij{1,t}=[a_valij{1,t},c_valij{1,t},1];
    dynamic_blx{1,t}=[lower_x{1,t};lower_u{1,t};-inf];
    dynamic_bux{1,t}=[upper_x{1,t};upper_u{1,t};inf];
    dynamic_blc{1,t}=[zeros(size_b(t),1);-big_m];
    dynamic_buc{1,t}=[zeros(size_b(t),1);inf];
    dynamic_c{1,t}=[costc{1,t};costd{1,t};1];
end

dynamic_subi{1,T}=[a_subi{1,T},c_subi{1,T},];
dynamic_subj{1,T}=[a_subj{1,T},c_subj{1,T}+size_x(T)];
dynamic_valij{1,T}=[a_valij{1,T},c_valij{1,T}];
dynamic_blx{1,T}=[lower_x{1,T};lower_u{1,T}];
dynamic_bux{1,T}=[upper_x{1,T};upper_u{1,T}];
dynamic_blc{1,T}=[zeros(size_b(T),1)];
dynamic_buc{1,T}=[zeros(size_b(T),1)];
dynamic_c{1,T}=[costc{1,T};costd{1,T}];

Cum_Probas=cell(1,T-1);
for t=1:T-1
    Cum_Probas{1,t}=[0,cumsum(probabilities{1,t})];
end

End_Algo=1;
iter=1;
Costs=[];

while End_Algo
    iter
    tic
    total_cost=0;
    trial_states=cell(1,T-1);
    clear prob;
    prob.c=dynamic_c{1,1};
    prob.blx=dynamic_blx{1,1};
    prob.bux=dynamic_bux{1,1};
    prob.buc=dynamic_buc{1,1};
    prob.blc=dynamic_blc{1,1};
    prob.a=sparse(dynamic_subi{1,1},dynamic_subj{1,1},dynamic_valij{1,1},size_b(1)+iter,size_x(1)+size_u(1)+1);
    [~,res]=mosekopt('minimize echo(0)',prob);
    sol=res.sol.bas.xx;
    solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
        disp('Unfeasible primal problem');
        pause
    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
        disp('Primal infinite optimal value');
        pause
    else
        trial_states{1,1}=sol(1:size_x(1));
    end
    zinf=sol'*dynamic_c{1,1};
    lower_bounds=[lower_bounds,zinf];
    total_cost=total_cost+dynamic_c{1,1}(1:size_x(1)+size_u(1))'*sol(1:size_x(1)+size_u(1));
    
    for t=2:T-1
        clear prob;
        prob.c=dynamic_c{1,t};
        prob.blx=dynamic_blx{1,t};
        prob.bux=dynamic_bux{1,t};
        Alea_Uniform=rand;
        [~,Index] = histc(Alea_Uniform,Cum_Probas{1,t-1});
        if (Alea_Uniform==1)
            Index=nb_scenarios_rhs(t);
        end
        aux1=bs{1,t}{1,Index}-sparse(b_subi{1,t-1},b_subj{1,t-1},b_valij{1,t-1},size_b(t),size_x(t-1))*trial_states{1,t-1};
        dynamic_buc{1,t}(1:size_b(t))=[aux1];
        dynamic_blc{1,t}(1:size_b(t))=[aux1];
        prob.buc=dynamic_buc{1,t};
        prob.blc=dynamic_blc{1,t};
        prob.a=sparse(dynamic_subi{1,t},dynamic_subj{1,t},dynamic_valij{1,t},size_b(t)+iter,size_x(t)+size_u(t)+1);
        [~,res]=mosekopt('minimize echo(0)',prob);
        sol=res.sol.bas.xx;
        solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
        if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
            disp('Unfeasible primal problem');
            pause
        elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
            disp('Primal infinite optimal value');
            pause
        else
            trial_states{1,t}=sol(1:size_x(t));
            total_cost=total_cost+dynamic_c{1,t}(1:size_x(t)+size_u(t))'*sol(1:size_x(t)+size_u(t));
        end
    end
    
    clear prob;
    prob.c=dynamic_c{1,T};
    prob.blx=dynamic_blx{1,T};
    prob.bux=dynamic_bux{1,T};
    Alea_Uniform=rand;
    [~,Index] = histc(Alea_Uniform,Cum_Probas{1,T-1});
    if (Alea_Uniform==1)
        Index=nb_scenarios_rhs(T);
    end
    aux1=bs{1,T}{1,Index}-sparse(b_subi{1,T-1},b_subj{1,T-1},b_valij{1,T-1},size_b(T),size_x(T-1))*trial_states{1,T-1};
    dynamic_buc{1,T}(1:size_b(T))=[aux1];
    dynamic_blc{1,T}(1:size_b(T))=[aux1];
    prob.buc=dynamic_buc{1,T};
    prob.blc=dynamic_blc{1,T};
    prob.a=sparse(dynamic_subi{1,T},dynamic_subj{1,T},dynamic_valij{1,T},size_b(T),size_x(T)+size_u(T));
    [~,res]=mosekopt('minimize echo(0)',prob);
    sol=res.sol.bas.xx;
    solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
        disp('Unfeasible primal problem')
        pause
    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
        disp('Primal infinite optimal value');
        pause
    else
        total_cost=total_cost+dynamic_c{1,T}(1:size_x(T)+size_u(T))'*sol(1:size_x(T)+size_u(T));
    end
    
    Costs=[Costs;total_cost];
    zsup=mean(Costs)+talpha*sqrt(var(Costs))/sqrt(iter);
    upper_bounds=[upper_bounds,zsup];
    
    %Backward pass
    %t=2 slope beta_2^k intercept theta_2^k
    
    for t=T:-1:2
        if (t==T)
            intercept=0;
            slope=zeros(size_x(T-1),1);
            for j=1:nb_scenarios_rhs(T)
                clear prob;
                aux1=bs{1,T}{1,j}-(sparse(b_subi{1,T-1},b_subj{1,T-1},b_valij{1,T-1},size_b(T),size_x(T-1)))*trial_states{1,T-1};
                prob.buc=[aux1];
                prob.blc=[aux1];
                prob.c=dynamic_c{1,T};
                prob.blx=dynamic_blx{1,T};
                prob.bux=dynamic_bux{1,T};
                prob.a = sparse(dynamic_subi{1,T},dynamic_subj{1,T},dynamic_valij{1,T},size_b(T),size_x(T)+size_u(T));
                [~,res]=mosekopt('minimize echo(0)',prob);
                sol=res.sol.bas.xx;
                solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
                if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                    disp('Unfeasible primal problem')
                    pause
                elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                    disp('Primal infinite optimal value');
                    pause
                else
                    dual1=res.sol.bas.slc;
                    dual2=res.sol.bas.suc;
                    dual3=res.sol.bas.slx;
                    dual4=res.sol.bas.sux;
                    auxs=-sparse(b_subj{1,T-1},b_subi{1,T-1},b_valij{1,T-1},size_x(T-1),size_b(T))*(dual1(1:size_b(T))-dual2(1:size_b(T)));
                    aux=(dual1(1:size_b(T))-dual2(1:size_b(T)))'*prob.buc(1:size_b(T))-auxs'*trial_states{1,T-1};
                    for i=1:(size_x(T)+size_u(T))
                        if (prob.blx(i)~=-inf)
                            aux=aux+prob.blx(i)*dual3(i);
                        end
                        if (prob.bux(i)~=inf)
                            aux=aux-prob.bux(i)*dual4(i);
                        end
                    end
                    slope=slope+probabilities{1,T-1}(j)*auxs;
                    intercept=intercept+probabilities{1,T-1}(j)*aux;
                end
            end
        else
            intercept=0;
            slope=zeros(size_x(t-1),1);
            for j=1:nb_scenarios_rhs(t)
                clear prob;
                prob.c=dynamic_c{1,t};
                prob.blx=dynamic_blx{1,t};
                prob.bux=dynamic_bux{1,t};
                aux1=bs{1,t}{1,j}-sparse(b_subi{1,t-1},b_subj{1,t-1},b_valij{1,t-1},size_b(t),size_x(t-1))*trial_states{1,t-1};
                dynamic_buc{1,t}(1:size_b(t))=[aux1];
                dynamic_blc{1,t}(1:size_b(t))=[aux1];
                prob.buc=dynamic_buc{1,t};
                prob.blc=dynamic_blc{1,t};
                prob.a=sparse(dynamic_subi{1,t},dynamic_subj{1,t},dynamic_valij{1,t},size_b(t)+iter+1,size_x(t)+size_u(t)+1);
                [~,res]=mosekopt('minimize echo(0)',prob);
                sol=res.sol.bas.xx;
                solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
                if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                    disp('Unfeasible primal problem');
                    pause
                elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                    disp('Primal infinite optimal value');
                    pause
                else
                    dual1=res.sol.bas.slc;
                    dual2=res.sol.bas.suc;
                    dual3=res.sol.bas.slx;
                    dual4=res.sol.bas.sux;
                    auxs=-sparse(b_subj{1,t-1},b_subi{1,t-1},b_valij{1,t-1},size_x(t-1),size_b(t))*(dual1(1:size_b(t))-dual2(1:size_b(t)));
                    aux=(dual1(1:size_b(t))-dual2(1:size_b(t)))'*prob.buc(1:size_b(t))+ dual1(size_b(t)+1:size_b(t)+iter+1)'*prob.blc(size_b(t)+1:size_b(t)+iter+1)-auxs'*trial_states{1,t-1};
                    for i=1:(size_x(t)+size_u(t))
                        if (prob.blx(i)~=-inf)
                            aux=aux+prob.blx(i)*dual3(i);
                        end
                        if (prob.bux(i)~=inf)
                            aux=aux-prob.bux(i)*dual4(i);
                        end
                    end
                    slope=slope+probabilities{1,t-1}(j)*auxs;
                    intercept=intercept+probabilities{1,t-1}(j)*aux;
                end
            end
        end
        dynamic_subi{1,t-1}=[dynamic_subi{1,t-1},(size_b(t-1)+iter+1)*ones(1,size_x(t-1)+1)];
        dynamic_subj{1,t-1}=[dynamic_subj{1,t-1},[1:1:size_x(t-1)],size_x(t-1)+size_u(t-1)+1];
        dynamic_valij{1,t-1}=[dynamic_valij{1,t-1},-slope',1];
        dynamic_buc{1,t-1}=[dynamic_buc{1,t-1};inf];
        dynamic_blc{1,t-1}=[dynamic_blc{1,t-1};intercept];
    end
    time=[time;toc];
    End_Algo=abs((zsup-zinf)/zsup)>tol;
    if (iter>=nb_iter_max)
        End_Algo=0;
    end
    iter=iter+1;
end



