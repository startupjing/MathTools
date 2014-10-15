classdef simplexMethod < handle
    % The Matlab code implements a class that can solve linear programming
    % problems using big-M method
    % 
    % The problem is in the following format
    % max(or min) c'x
    % subject to Ax>=b or Ax<=b or Ax=b
    %            x>=0
    % Note that RHS b should be non-negative
    % If it is negative, convert it by multiplying -1 on both sides
    % --------------------------
    % Example usage
    %    max 2x1 ? 6x2
    % s.t.  ?x1 ? x2 ?x3 ? ?2
    %       2x1 ? x2 + x3 ? 1 
    %       x1 , x2 , x3 ? 0
    % Since RHS of ?x1 ? x2 ?x3 ? ?2 is -2 and is negtive,
    % we conver it to x1 + x2 + x3 >= 2
    % 
    % Now the problem becomes
    %    max 2x1 ? 6x2
    % s.t.  x1 + x2 + x3 ? 2
    %       2x1 ? x2 + x3 ? 1 
    %       x1 , x2 , x3 ? 0
    % 
    % We should define c = [2 -6 0]' , A = [1 1 1; 2 -1 1] 
    % b = [2 1]' , isEq = [1 -1]
    % Then we construct an instance by
    % solver = simplexMethod(A,b,c,isEq,'max');
    % and solver.compute will give you the result
    % ---------------------------
    
    % user-defined properties
    properties
        A;     % coefficient matrix for constraints
        b;     % RHS matrix
        c;     % objective coefficients
        isEq;  % vector that identifies <= (use -1) or >= (use 1)
        m;     % number of constraints
        n;     % number of variables
        opt;   % final optimal solution
        val;   % final optimal value
        option;  % choose 'max' or 'min'
        M;    % big M
        slack_num;  % number of slack variables
        times;  % current number of iterations
        z_cof;  % coefficients in objective function
        stop;   % terminate flag
        var;    % variables
        num_art;    % number of artificial variables
        remove_art;  % removed artificial variables
        art_idx;     % index of artificial variables
        basic;      % basis variables
        select_row;  % seleted row number
        tab_form;    % tabular matrix
        var_unbound;  % unbounded variable
        choice;  % 1 for optimal, 2 for unbounded, 3 for infeasible, 4 for max_iteration exceeds
    end
    
    
    
    methods
        % constructor
        function p = simplexMethod(A,b,c,isEq,option)
            % check number of input
            if nargin < 4
                error('Number of inputs is incorrect');
            end
            
            % check option
            % change the sign of coefficients in objective function
            if strcmpi(option,'max')
                c = -c;
            end
            
            % check if RHS is negative
            if any(b) < 0
                error('RHS canot be negative. Please convert to correct format');
            end
            
            % set fields
            p.A = A;
            p.b = b;
            p.c = c;
            p.times = 0;
            p.isEq = isEq;
            p.option = option;
            p.stop = 0;
            % m is number of constraints 
            % n is number of variables
            [p.m, p.n] = size(A);
            
            % set value for bigM
            p.M = ceil(max(norm(A),norm(b))/100)*100;
            
            % build variable list {x1,x2, ... xn}
            p.var = cell(1,p.n);
            % build strings for variables
            % in the form of "x1, x2 ..."
            for i = 1:p.n
                p.var{i} = ['x' num2str(i)];
            end   
        end
        
        % actual function that computes optimal solution
        % will give information if infeasible or unbounded
        function compute(p)
            p.addSlack;
            p.prepare;
            p.printTab;
            while ~p.stop
                p.iterate;
                p.times = p.times + 1;
                p.printTab;
                p.check_art;
                p.is_end;
            end;
            p.printSoln;
        end
        
        % add slack variable
        function addSlack(p)
            % count number of slack variables
            % need to add slack variables if it is inequality
            p.slack_num = sum(p.isEq < 0) + sum(p.isEq > 0);
            % enlarge matrix A and objective row c
            p.A = [p.A zeros(p.m, p.slack_num)];
            p.c = [p.c; zeros(p.slack_num,1)];
            % starting index of slack variables
            id = p.n + 1;
            % loop every row, set coefficients of slack variables
            for i = 1:p.m
                % set -1 if LHS>=RHS and set 1 if LHS<=RHS
                p.A(i, id) = -p.isEq(i);
                % go to next slack variable
                id = id + 1;
                % build string for slack varible
                % in the form of "s1, s2 ..."
                p.var{end+1} = ['s' num2str(i)];
            end    
        end
        
        % build the initial tabular form
        function prepare(p)
            % loop through every variable
            % find basis variables (slack variable with coefficient 1)
            for i = 1 : p.n + p.slack_num
                % check if the variable column only contains
                % one non-zero element
                if nnz(p.A(:,i)) == 1
                    % find row number of value 1
                    % if the variablecolumn is 0..0 1 0 ...0
                    row = find(p.A(:,i) == 1);
                    % if the only non-zero element in the variable column
                    % is 1, do the following
                    if ~isempty(row)
                        % add row number to seleted_row list 
                        % if it is not in the list previously
                        if ~ismember(row,p.select_row)
                            p.select_row(end+1) = row;
                            % add the variable to basis variable list
                            % i means the ith variable
                            p.basic(row) = i;
                        end
                    end
                end
            end
            
            % find number of artificial variable needed
            p.num_art = p.m - sum(p.basic > 0);
            % add artificial variable if needed
            if p.num_art > 0
                % set remove flag
                p.remove_art = 0;
                % index list for artificial variables
                p.art_idx = (p.n+p.slack_num+1):(p.n+p.slack_num+p.num_art);
                % add blocks for artificial variables
                p.A = [p.A zeros(p.m, p.num_art)];
                % find the row number where artificial variable is needed
                % in other words, the corresponding initial contraint is in
                % the form LHS - slack = RHS(positive)
                art_row = setdiff(1:p.m, p.select_row);
                
                
                % insert artificial variables
                for i = 1:length(art_row)
                    % enlarge matrix A to add columns for artificial
                    % variables and set its coefficient to be 1
                    p.A(art_row(i), p.n + p.slack_num + i) = 1;
                    % build string for artificial variables
                    % in the form "a1 a2 ... "
                    p.var{end+1} = ['art' num2str(i)];
                    % add artificial variable to basis
                    p.basic(art_row(i)) = p.n + p.slack_num + i;
                end
                % rebuild objective row
                p.c = [p.c; p.M*ones(p.num_art,1)];
            end
            
            % initialize coefficients in objective function
            p.z_cof = zeros(1, p.n + p.slack_num + p.num_art);
            % get coefficients for basis
            y = p.c(p.basic);
         
            % construct objective function coefficients
            p.z_cof = (p.c - p.A'*y)';
          
            % build tabular form
            p.tab_form = [p.A p.b];
            p.tab_form = [p.tab_form; [ p.z_cof -p.b'*y ]];
        end
        
        % perform one iteration of simplex method
        function iterate(p)
            % find column number of
            % smallest cofficient in objective row
            [~, col] = min(p.tab_form(end,1:end-1));
            
            % find row number of positive elements
            % in the previously selected column
            pos = find(p.tab_form(1:end-1, col)>0);
            
            % find index of smallest positive ratio
            [~, idx_ratio] = min(p.tab_form(pos,end)./p.tab_form(pos,col));
            idx_ratio = pos(idx_ratio);
            
            % update basis
            p.basic(idx_ratio) = col;
            
            % find pivot element
            pivot = p.tab_form(idx_ratio, col);
            
            % use division to set pivot to 1
            p.tab_form(idx_ratio, :) = p.tab_form(idx_ratio,:)/pivot;
            
            % find index of rows whose corresonding element is not 1
            row_idx = setdiff(1:p.m+1, idx_ratio);
            
            % update the table by row reduction
            p.tab_form(row_idx,:) =... 
               p.tab_form(row_idx,:) -...
               ones(p.m,1)*p.tab_form(idx_ratio,:) .*...
               (p.tab_form(row_idx,col)*ones(1,size(p.tab_form,2)));
        end
        
        % remove artificial variable if necessary
        function check_art(p)
            % if artificial variables have not been romoved
            % do the following
            if ~p.remove_art
                % remove artificial variable if its coefficient is positive
                if all(p.tab_form(end,p.art_idx) > 0)
                    % clear artificial variables in table
                    p.tab_form(:,p.art_idx) = [];
                    % clear them in variable list
                    p.var(p.art_idx) = [];
                    % set remove flag
                    p.remove_art = 1;
                    p.printTab;
                end
            end
        end
        
       
        % check if simplex method can be terminated
        % and decide optimal/infeasible/unbounded solution
        function is_end(p)
            % set terminate flag
            p.stop = 0;
            
            % terminate the algorithm if using too many iterations
            % (exceeds the maximum allowed number of iterations)
            if p.times >= 30
                p.stop = 1;
                % choice 4 is exceeding max_iter
                p.choice = 4;
            end   
            
            % terminate if all coefficients in the objective row
            % are non-negative
            if all(p.tab_form(end, 1:end-1) > -1e-10)
                % choice 1 is optimal
                p.stop = 1;
                
                % for big-M method, solution is infeasible
                % if one of artificial variable is non-zero
                % in other words, it has not been removed 
                % by function check_art(p)
                if ~p.remove_art
                    % choice 3 is infeasible
                    p.choice = 3;
                else
                    p.choice = 1;
                    p.opt = zeros(p.n+p.slack_num+p.num_art,1);
                    p.opt(p.basic) = p.tab_form(1:end-1,end);
                    if strcmpi(p.option, 'max')
                        p.val = p.tab_form(end,end);
                    elseif strcmpi(p.option, 'min')
                        p.val = -p.tab_form(end,end);
                    end
                end
         
            else
                % find column number in the objective row
                % where the coefficient is negative
                neg_col = find(p.tab_form(end,1:end-1)<0);
                % find number of negative elements in that column
                idx_unb = sum(p.tab_form(1:end-1, neg_col)<0);
                % check if all elements in that column are negative
                % in other words, check to see if there is no leaving
                % variable
                idx_unb = (idx_unb==p.m);
                % for big-M method, solution is unbounded if no leaving
                % variable exists
                if any(idx_unb)
                    % choice 2 is unbounded
                    p.choice = 2;
                    p.stop = 1;
                    p.var_unbound = neg_col(idx_unb>0);
                end
            end
        end
            
        % print tabular form for each iteration
        function printTab(p)
            fprintf('\n\n');
            % print iteration info
            fprintf('========== Iteration %2s ==========\n\n', num2str(p.times));

            % Print variables information
            for k=1:length(p.var)+1
                if k==1
                    fprintf('BV |', '');
                else
                    fprintf('%7s %4s', p.var{k-1}, '|');
                    if k == length(p.var)+1
                        fprintf(' %7s', 'RHS');
                    end
                end
            end
            fprintf('\n');
            fprintf('------------------------------------------------------------------\n');
            
            % print objective row
            fprintf('%2s | ', 'Z');
            for j = 1:size(p.tab_form,2)
               fprintf('%7.2f %4s', p.tab_form(p.m+1,j), '');
               if j == size(p.tab_form,2)-1
                       fprintf(' | ') 
               end
            end
            fprintf('\n');
            
            
            % print other rows
            for i=1:p.m
                fprintf('%2s | ', p.var{p.basic(i)});
                for j = 1:size(p.tab_form,2)
                    fprintf('%7.2f %4s', p.tab_form(i,j), '');
                    if j == size(p.tab_form,2)-1
                       fprintf(' | ') 
                    end
                end
                fprintf('\n');
            end
            
        end
        
        % print final result
        function printSoln(p)
            fprintf('\n\n');
            if p.stop
                fprintf('======== Result ========. \n');
                switch lower( p.choice )
                    % optimal solution
                    case 1
                        fprintf('Solution is optimal: \n');
                        for i=1:p.n
                            fprintf('x[%2d ] = %5.2f\n', i, p.opt(i));
                        end
                        fprintf('Optimal value: %5.2f\n', p.val);
                    % unbounded case
                    case 2
                        fprintf('Unbouned problem. \n');
                    % infeasible case
                    case 3
                        fprintf('Infeasible problem.\n');
                    % too many iterations
                    case 4
                        fprintf('Too many iterations.\n');
                        fprintf('Please Check input or algorithm.\n');
                end
                
            end
        end
    end
                
                
                    
        
        
end

