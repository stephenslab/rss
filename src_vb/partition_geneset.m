function sumstat_cell = partition_geneset(sumstat)
% USAGE: partition the summary data in a gene set to cell arrays by chromosome
% INPUT:
%	sumstat: a structure containing summary data (stored as array)
% OUTPUT:
%	sumstat_cell: a structure containing summary data (stored as cell array)
 
  % load the summary data stored as array
  betahat_mat = sumstat.betahat;
  se_mat      = sumstat.se;
  SiRiS_mat   = sumstat.SiRiS;
  chr_mat     = sumstat.chr;
  pos_mat     = sumstat.pos;

  % create the cell arrays
  chr_set  = unique(chr_mat);
  num_cell = length(chr_set);
  betahat  = cell(num_cell, 1);
  se       = cell(num_cell, 1);
  SiRiS    = cell(num_cell, 1);
  chr      = cell(num_cell, 1);
  pos      = cell(num_cell, 1);

  % make sure the chr order is ascending
  if ~issorted(chr_set)
    error('Error in chromosome order.')
  end

  % assign data to cell arrays
  for c=1:num_cell
    tmp_index = (chr_mat == chr_set(c));
 
    betahat{c,1} = betahat_mat(tmp_index, 1);
    se{c,1}      = se_mat(tmp_index, 1);
    SiRiS{c,1}   = SiRiS_mat(tmp_index, tmp_index);
    chr{c,1}     = chr_mat(tmp_index, 1);
    pos{c,1}     = pos_mat(tmp_index, 1);

    if ~issorted(pos{c,1})
      error('Error in base pair.')
    end 
  end

  % sanity check
  check_betahat = all(betahat_mat == cell2mat(betahat));
  check_se      = all(se_mat == cell2mat(se));
  check_chr     = all(chr_mat == cell2mat(chr));
  check_pos     = all(pos_mat == cell2mat(pos));
  check_all     = all([check_betahat check_se check_chr check_pos]);
  if ~check_all
    error('Error in partitioning the summary data of the gene set.')
  end

  % save the cell arrays in a structure
  sumstat_cell.betahat = betahat;
  sumstat_cell.se      = se;
  sumstat_cell.SiRiS   = SiRiS;
  sumstat_cell.chr     = chr;
  sumstat_cell.pos     = pos;

end 
