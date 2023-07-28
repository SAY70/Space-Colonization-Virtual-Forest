classdef SC_tree_obj < matlab.mixin.Copyable % shallow copy [can use copy()]
    
    properties
    
        tbuds               % terminal buds (nodes
        b_m     % = 0.5  ;  % branch length magnitude (for each node segment)
        min_d   % =  1 ;    % minimum distance to branch
        max_d   % =  3 ;    % maximum distance to branch
        max_its % = 20 ;    % max num of iterations before quitting fx
        trunk               % trunk information
        tree                % tree information (trunk + branches)
        max_num_brnch = 5e3;% maximum number of branches for tree
        
        % scobi parameters
        pos             % tree position ((head-tail)/2)
        diel_r          % real of complex dielectric [eps = eps' + i(eps'')]    [F/m]
        diel_i          % imag of complex dielectric [eps = eps' + i(eps'')]    [F/m]
        semi_minor      % semi-minor axis of cylinder (branch) or disc (leaf)   [m]
        semi_major      % semi-major axis of cylinder (branch) or disc (leaf)   [m]
        len             % length of cylinder or disc                            [m]
        type string = '' % type (trunk, branch, leaf)
        
    end
    
    methods
        
        function set_inputs(obj,tbuds,b_m,min_d,max_d,max_its)
            
            obj.tbuds = tbuds;
            obj.b_m = b_m;
            obj.min_d = min_d;
            obj.max_d = max_d;
            obj.max_its = max_its;
            
        end
        
        function [heads, tails] = make_root(obj)

            get_u = @(a) a./vecnorm(a);

            heads = []   ;
            tails = []   ;

            c_tail = [0 0 0]';
            c_head = [0 0 1]' * obj.b_m;
            c_uvec = get_u(c_head-c_tail); 

            found = false;

            heads = [heads c_head] ;
            tails = [tails c_tail] ;
            
            its = 0;

            % make trunk    
            while ~found && its < obj.max_its

                % current distance / magnitude
                cdm = vecnorm(obj.tbuds - c_head);
                % closest leaves
                inds = cdm<obj.max_d ;

                if ~any(inds)

                    c_tail = c_head;
                    c_head = c_head + c_uvec * obj.b_m;

                    % add to matrix
                    heads = [heads c_head] ;
                    tails = [tails c_tail] ;
                    its = its+1;
                else
                    % trunk is in close proximity to the buds; escape fx
                    found = true ;
                    its = 0;
                end  
            end
            
            obj.trunk.heads = heads ;
            obj.trunk.tails = tails ;


        end
        
        function [bh,bt] = make_branch(obj)

            % bh -- branch heads
            % bt -- branch tails
            % tbuds -- terminal buds of tree branches
            % l -- lower bound / least distance
            % u -- upper bound / maximum distance
            % b_m -- branch magnitude
            % max_its -- maximum number of iterations allowed where we grow but
            %            do not find a bud

            % for each terminal bud, determine if there is a branch
            % that is between the least and maximum distance allowable

            get_u = @(a) a./vecnorm(a);
            plt3r = @(a) plot3(a(1,:), a(2,:), a(3,:), 'color', 'r');
        %     count = zeros([1 size(tbuds,2)]);
        
            bh = obj.trunk.heads ;
            bt = obj.trunk.tails ;

            dbuds  = obj.tbuds;                       % deletable bud array
            cl_br  = zeros( [1, size(obj.tbuds,2)] ); % closest branch
            dists  = zeros( [1, size(obj.tbuds,2)] ); % distances
            uvecs  = zeros( [3, size(obj.tbuds,2)] ); % unit vectors 

            its = 0;
        %     max_its = 20;

            % turn this into a fx?
            while ~isempty(dbuds) && (its < obj.max_its) && (size(bh,2)<obj.max_num_brnch);

            for i = 1:size(dbuds,2)

                % current bud
                cur_bud = dbuds(:,i) ;
                % all distances of each branch to current bud
                dist = vecnorm( cur_bud - bh ) ;
                md = min(dist) ;
                mdi = find( dist == min(dist) );

                dists(i) = md  ;
                mdi = mdi(1) ; % in the event a point has identical distances
                cl_br(i) = mdi ;
                uvecs(:,i) = get_u( cur_bud - bh(:, mdi) );


            end

            % unique distances
            ucl_br = unique(cl_br) ;
            nh = zeros([3, numel(ucl_br)]);
            nt = nh;
            ii = 0;
            for i =1:length(ucl_br)
                % branch indexes
                br_i = ucl_br(i) ;
                % unit vector indexes
                ci = find( cl_br == ucl_br(i) ) ;
                % get average unit vector
                auvec = mean( uvecs(:, ci), 2 ) ;
                % new unit vector, averaged by counts
                nv = auvec * obj.b_m ;
                new_tail = bh(:,br_i) ;
                nv = new_tail + nv ;

                % update index
                ii = ii + 1;
                nt(:,ii) = new_tail ;
                nh(:,ii) = nv ;

        %         plt3r([ new_tail nv ])

            end

            % add new branches to total list (cannot preallocate)
            bh = [bh nh] ;
            bt = [bt nt];

            % delete any leaf that has been reached
            if any(dists < obj.min_d)
                dbuds = obj.del_a3( dbuds, find(dists < obj.min_d) );
                % reset kill count
                its = 0 ;
            else 
                % increment kill count
                its = its + 1 ;

            end
            end
            
            obj.tree.heads = bh;
            obj.tree.tails = bt;



        end
        
        function shinozaki_pipe(obj, dbh, d0)
            
            % an implementation based on the Shinozaki pipe model theory of
            % plant cylinder construction.
            % https://doi.org/10.18960/seitai.14.3_97
            % with implementation based on the description of 
            % https://algorithmicbotany.org/papers/colonization.egwnp2007.pdf
            %
            % instead of using an adjustment parameter [Runions 2007, eq.1],
            % pipes are formed by determining the total number of children 
            % for the plant. the counts are then scaled to fit between a
            % maximum of the diameter at breast height and the minimum
            % desired radius for terminal branch nodes.
            %
            % this code is designed for cylindrical branches with equal
            % semi-major and semi-minor axes. space is made to compute
            % semi-major and semi-minor independently, though i'm sure
            % there are several considerations for the orientation of the
            % cylinder along the growing structure.
            %
            % inputs %
            % dbh -- diameter at breast height [m]
            % d0  -- diameter of the terminal branches [m]
            %% find terminal branches


            h = obj.tree.heads;
            t = obj.tree.tails;

            % lambda functions
            ufndc = @(a,b) find( sum( a ==  b, 1) == 3 );
            get_par   = @(ch) ufndc( ch, h ) ;
            get_child = @(ch) ufndc( ch, t ) ;

            ts = zeros([1 size(h,2)]); % terminal points
            % d0 = 0.5;             % radius of terminal branches [m]
            as = zeros(size(ts)); % semi-major axis [m] (typically a==b)
            bs = zeros(size(ts)); % semi-minor axis [m] (typically a==b)
            is = [];

            % find terminal branches
            for i = 1:length(h)
                if ~ismember(h(:,i)', t', 'rows')
                    is = [is i];
                end
            end

            % initialize terminal points with terminal diameter
            as(is) = d0;
            bs(is) = d0;

            next_cis = is ;
            max_its = 1000;
            its = 0;
            
            % determine diameter of all branches
            while any(as==0)
                % reset parameters on while loop
                cis = next_cis;
                next_cis = [] ;
                % loop over current end points
                for i =1:length(cis)
                    
                    % check the radius of current branch's children
                    cind = cis(i);
                    ch = t(:, cis(i)) ;
                    child = get_child( ch );

                    % if any of the children are not defined, then
                    % come back later
                    if any(as(child) == 0)
                        next_cis = [next_cis cind] ;
                    else
                    % if all children have defined diameter, calculate next
                    % diameter
                    par = get_par( ch );

                    new_a = 0;
                    new_b = 0;
                    for ii = 1:length(child)
                        cur_ci = child(ii) ;
                        new_a = new_a + as(cur_ci) ;
                        new_b = new_b + bs(cur_ci) ;
                    end

                    as( par ) = new_a ;
                    bs( par ) = new_b ;
                    
                    % this iteration's tail = next iteration's head
                    next_cis = [next_cis par] ;
                    end

                end
                
                % emergency break code
                % since we're in a while loop
                if isequal(cis, next_cis)
                    its = its + 1;
                end
                % break if we exceed desired iterations
                if its == max_its
                    warning('shinozaki pipe model ended unexpectedly')
                    break
                end

            end


            obj.semi_major = rescale( as, d0, dbh ) ;
            obj.semi_minor = rescale( bs, d0, dbh ) ;
        end
        
        function o = del_a3(~, A, i, rc)
            % delete column/row vector from array

            % A -- input array (2d)
            % i -- index to delete
            % rc -- row or column?
            %       1 == column
            %       2 == row   
            % o -- output array

            % ex
            % A = reshape( 1:3*12, [3,12]);
            % i = [4 9];

            if ~exist('rc', 'var')
                rc = 1 ; % column
            end

            % unfortunate, size-changing for loop but it's for the best
            ii = [] ;
            for n = 1:size(A,2)
                if ~any( ismember(n,i) )
                    ii = [ii n];
                end
            end

            % delete element from array
            switch rc
                case 1
                    o = A(:,ii) ;
                case 2
                    o = A(ii,:) ; 
            end

        end
        
        function assign_diel(obj)
            
            trunk_inds  = (obj.type == "trunk" ) ;
            branch_inds = (obj.type == "branch") ;
            
            obj.diel_r(trunk_inds)  = 15.6 ; % [F/m]
            obj.diel_i(trunk_inds)  =  3.8 ; % (+ sqrt(i)) 
            obj.diel_r(branch_inds) = 12.0 ; % [F/m] 
            obj.diel_i(branch_inds) = 2.93 ; % (+ sqrt(i)) 
            
        end
        
        function vol = get_vol(obj)
            % volume of cylinders
            vol = pi .* obj.semi_major .* ...
                obj.semi_minor .* obj.len ; 
        end
        
        function get_pos(obj)
            % return position [m]
            obj.pos = (1/2) .* ( obj.tree.heads - obj.tree.tails) ;
        end
        
        function consolidate_vectors(obj)
            % L-systems will usually make a vector "long" by stacking
            % smaller vectors together. This method is a best effort
            % attempt at taking several vectors that are connected from
            % head to tail in the same unit vector direction and
            % consolidating them into a single vector with a new magnitude
            %
            % please only use this function when vectors are of the same
            % magnitude.
            
            % lambda functions
            % fndc  = @(a,b) find( sum( a == (a - b), 1) == 3 );
            ufndc = @(a,b) find( sum( a ==  b, 1) == 3 );
            nmlz  = @(a) round( a ./ vecnorm(a), 10);
            rvn   = @(a) round( vecnorm(a), 10) ;
            
            % TO-DO: preallocat with nans and filter out instead of
            % constant array resizing
            n_tails = [] ;
            n_heads = [] ;

            mag = obj.b_m ;
            heads = obj.tree.heads ;
            tails = obj.tree.tails ;
            
            for c_i1 = 1 : size(heads,2)
                
                c_h1 = heads(:,c_i1);
                c_t1 = tails(:,c_i1);
                % current unit vector
                uvec1 = nmlz( c_h1 - c_t1 );
                % all matching unit vectors in list
                uvec2 = nmlz( heads - c_t1 );
                % get index of all matching unit vectors
                c_i2 = ufndc(uvec1,uvec2);

                % get subset of all head and tail vectors for testing
                theads = heads(:,c_i2) ;
                ttails = tails(:,c_i2) ;

                % the current head and tail
                ct_tail = c_t1;
                ct_head = c_h1;
                ct_tail2 = c_t1;
                
                % iterate over the vectors being tested
                for ti = 1:size(theads,2)
                    tmpi = find( rvn(theads-ct_tail) == ti*mag) ;
        %             disp(tmpi)
                    if ~isempty(tmpi)
                        if isequal(rvn(ttails(:,tmpi)-round(ct_tail2,14)), mag);
                            ct_tail2 = ttails(:,tmpi) ;
                            ct_head = theads(:,tmpi);
                        end
                    end
                end
                
                
                % assign values to output
                
                if isempty(n_heads)
                    % first iteration
                    n_tails = [ct_tail] ;
                    n_heads = [ct_head] ;
                elseif ~ismember(n_heads',ct_head','rows')
                    % only assign if the head vector is not present in the
                    % file. this could potentially cause a problem for two 
                    % branches terminating at the same point.
                    n_tails = [n_tails ct_tail] ;
                    n_heads = [n_heads ct_head] ;
                end
                
                
            end
            
            % assign updated values to object
%             obj.tree.heads2 = n_heads;
%             obj.tree.tails2 = n_tails;
            
            obj.tree.heads = n_heads;
            obj.tree.tails = n_tails;
        end
         
        function identify_branches(obj)
            
            % update the vector length
            obj.len = vecnorm( obj.tree.heads - obj.tree.tails ) ;
            
            % identify vectors
            uvecs = ( obj.tree.heads - obj.tree.tails ) ./ obj.len ;
            logic1 = ( sum(uvecs == [0 0 1]',1) == 3 ) ;
            logic2 = ( sum(obj.tree.tails == [0 0 0]',1) == 3 ) ;
            logic = logic1 & logic2;
            obj.type(logic)  = 'trunk';
            obj.type(~logic) = 'branch';
            
            
        end
    end
    
end