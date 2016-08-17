function [rank_new, gamma_logratio] = propose_gamma(rank_old, p_gamma, repeat)
% USAGE: propose the inclusion indicator vector gamma based on the marginal rank
% SOURCE: based on Guan & Stephens (2011) Appendix A and Zhou et al (2013) Text S2
% INPUT:
%	rank_old: the rank of snps in the current model, 1 by ngamma_old
%	p_gamma: the pmf of the mixture of discrete uniform and geometric rv, p by 1
%       repeat: the number of local proposals to compound required by small world, scalar
% OUTPUT:
%       rank_new: the rank of snps included in the updated model, 1 by ngamma_new
%	gamma_logratio: log ratio contribution of gamma to the MH log ratio, scalar
	
	p = length(p_gamma); % the total number of snps

	gamma_logratio 	= 0;
	ngamma_old 	= length(rank_old); % the number of snps included in the current model
	ngamma_new 	= ngamma_old;
	rank_new 	= rank_old;
	
	% create a map object for snps in the current model
	snpmap = containers.Map('KeyType', 'double','ValueType', 'uint8');

	% initialize it if ngamma_old is not zero
	if ngamma_old ~= 0
		k 	= rank_old; 				% keys: the rank of snps, double
		v 	= ones(1, ngamma_old, 'uint8'); 	% mapped values: 1 (unit8)
		snpmap 	= containers.Map(k, v);  	
	end
	
	for i = 1:repeat
		
		% rare case 1: no snp in the current model; must add one snp
		if ngamma_new == 0
			r_add 		= sample(p_gamma);
			snpmap(r_add) 	= 1;
			rank_new 	= r_add;
			ngamma_new 	= ngamma_new + 1;		
			gamma_logratio 	= gamma_logratio - log(p_gamma(r_add)); % compute the contribution to log MH ratio
			continue 						% skip the rest of current iteration 
		end

		% rare case 2: all snps in the current model; must remove one snp
		if ngamma_new == p 
			col_id 		= randint(1,ngamma_new); 			% uniformly pick one snp to remove
			r_remove 	= rank_new(col_id);
			remove(snpmap, r_remove);
			rank_new(col_id) = []; 					% remove the snp	
			gamma_logratio 	 = gamma_logratio + log(ngamma_new); 	% compute the contribution to log MH ratio
			ngamma_new 	 = ngamma_new - 1;
			continue 						% skip the rest of current iteration
		end

		gamma_flag = rand; % decide the type of move
		
		if gamma_flag < 0.4 % add a snp
			% propose a snp NOT in the current model
			% randomly pick a snp until the picked snp is NOT in the current model
			while true 
				r_add = sample(p_gamma);
  				if (isKey(snpmap, r_add) == 0); break; end
			end
			% calculate Pr(picking a snp NOT in the current model)
			prob_total = 1;
			prob_total = prob_total - sum(p_gamma(rank_new));
			% add the proposed snp
			snpmap(r_add) 	= 1;
			rank_new 	= [rank_new r_add];
			ngamma_new 	= ngamma_new + 1;
			% calculate the contribution to log MH ratio
			gamma_logratio  = gamma_logratio - log(p_gamma(r_add) / prob_total) - log(ngamma_new);

		elseif gamma_flag >=0.4 && gamma_flag < 0.8 % remove a snp
			% uniformly draw a snp from the current model
			col_id 	 = randint(1,ngamma_new); 
			r_remove = rank_new(col_id); % return the location (a.k.a. rank) of the chosen snp
			% calculate Pr(picking a snp NOT in the proposed model)
			prob_total = 1;
			prob_total = prob_total - sum(p_gamma(rank_new));
			prob_total = prob_total + p_gamma(r_remove); % b/c this snp is NOT included in the proposed model
			% remove the proposed snp
			remove(snpmap, r_remove);
			rank_new(col_id) = [];
			% calculate the contribution to log MH ratio
			gamma_logratio 	= gamma_logratio + log(p_gamma(r_remove) / prob_total) + log(ngamma_new);
			ngamma_new 	= ngamma_new - 1;

		else % swap inclusion/exclusion status of two snps
			col_id 	 = randint(1,ngamma_new); % pick a snp in the model to remove
			r_remove = rank_new(col_id);
			while true
				r_add = sample(p_gamma); % randomly pick a snp based on marginal association rank
				if (isKey(snpmap, r_add) == 0); break; end % pick a snp NOT in the model to add
			end
			% calculate the probability of picking a snp NOT in the current model
			prob_total 	= 1;
			prob_total 	= prob_total - sum(p_gamma(rank_new));
			gamma_logratio 	= gamma_logratio + log( p_gamma(r_remove) / (prob_total + p_gamma(r_remove) - p_gamma(r_add) ) );
			gamma_logratio 	= gamma_logratio - log( p_gamma(r_add) / prob_total );
			% flip the indicators of these two snps
			remove(snpmap, r_remove);
			snpmap(r_add) 	 = 1;
			rank_new(col_id) = [];
			rank_new 	 = [rank_new r_add];
		end
	
	end
	
	% erase the map object to free up system memory
	clear('snpmap');
end
