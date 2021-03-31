library(MatrixEQTL)

### Format a number with comma delimited thousands
.s = function(x)formatC(x=x,digits=ceiling(log10(max(x)+1)),big.mark=',',big.interval=3);

### Single Simplified Matrix eQTL engine
.SingleMatrixEQTL <- setRefClass( ".SingleMatrixEQTL",
	fields = list( 
		zgene = "SlicedData",
		zsnps = "SlicedData",
		zcvrt = "matrix"
	),
	methods = list(
		initialize = function( snps, gene, cvrt ) {
			### Tests for consistency
			# Done before being called by the CorrMeta_testAll

			### Initialize
			{
				# Standard deviations before standardization
				gene.std = vector('list',gene$nSlices());
				snps.std = vector('list',snps$nSlices());
			}

			### Covariate Processing
			{
				# Combine and add constant element
				if( cvrt$nRows()>0 ) {
					cvrt$SetNanRowMean();
					cvrt$CombineInOneSlice();
					cvrt = rbind(matrix(1,1,snps$nCols()),cvrt$getSlice(1));
				} else {
					cvrt = matrix(1,1,snps$nCols());
				}

				# Standardize and orthogonolize via QR decomposition
				q = qr(t(cvrt));
				if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * snps$nCols() ) {
					stop("Colinear or zero covariates detected");
				}
				.self$zcvrt = t( qr.Q(q) );
			}
			
			### Expression
			{
				gene$SetNanRowMean();
				# Orthogonolize expression w.r.t. covariates
				for( sl in 1:gene$nSlices() ) { # sl = 1L
					slice = gene$getSlice(sl);
					rowsq1 = rowSums(slice^2);
					slice = slice - tcrossprod(slice,zcvrt) %*% zcvrt;
					rowsq2 = rowSums(slice^2);
					# kill rows colinear with the covariates
					delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
					if(any(delete.rows)) {
						slice[delete.rows,] = 0;
						rowsq2[delete.rows] = 1;
					}
					div = sqrt(rowsq2); #sqrt( rowSums(slice^2) );
		# 			div[ div == 0 ] = 1;
					gene.std[[sl]] = div;
					gene$setSliceRaw(sl, slice / div);
				}
				rm(rowsq1, rowsq2, delete.rows, div);
				rm( sl, slice );
				#gene$RowRemoveZeroEps();
			}
			
			### Genotype 
			{
				snps$SetNanRowMean();
				# Orthogonolize expression w.r.t. covariates
				# status("Orthogonolizing expression w.r.t. covariates");
				for( sl in 1:snps$nSlices() ) { # sl = 1L
					slice = snps$getSlice(sl);
					rowsq1 = rowSums(slice^2);
					slice = slice - tcrossprod(slice,zcvrt) %*% zcvrt;
					rowsq2 = rowSums(slice^2);
					# kill rows colinear with the covariates
					delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
					if(any(delete.rows)) {
						slice[delete.rows,] = 0;
						rowsq2[delete.rows] = 1;
					}
					div = sqrt(rowsq2); #sqrt( rowSums(slice^2) );
		# 			div[ div == 0 ] = 1;
					snps.std[[sl]] = div;
					snps$setSliceRaw(sl, slice / div);
				}
				rm(rowsq1, rowsq2, delete.rows, div);
				rm( sl, slice );
				#snps$RowRemoveZeroEps();
			}
			
			### Save preprocessed data sets
			.self$zsnps = snps;
			.self$zgene = gene;
		},
		get_tstats = function(geneslice) {

			dfFull = ncol(	zsnps ) - 1 - nrow(zcvrt);
			testfun = function(x) { return( x * sqrt( dfFull / (1 - pmin(x^2,1))));	}

			gslice = zgene[[geneslice]];
			
			rez = vector('list',zsnps$nSlices());
			for( sl in 1:zsnps$nSlices() ) { # sl = 1L
				rez[[sl]] = testfun(tcrossprod(gslice,zsnps[[sl]]));
			}
			return(rez);
		}
	)
)

.listBuilder1 = setRefClass(".listBuilder1",
	fields = list(
		dataEnv = "environment",
		n = "integer",
		count = 'numeric'
	),
	methods = list(
		initialize = function() {
			.self$dataEnv = new.env(hash = TRUE);
			.self$n = 0L;
			.self$count = 0;
# 			cumlength <<- 0;
			return(.self);
		},
		add = function(x) {
			if(length(x) > 0) {
				n <<- n + 1L;
# 				cumlength <<- cumlength + length(x);
				assign(paste(n), x, dataEnv );
				.self$count = .self$count + length(x);
			}
			return(.self);
		},
		set = function(i,x) {
			i = as.integer(i);
			if(length(x) > 0) {
				if(i>n)
					n <<- i;
				assign(paste(i), x, dataEnv );
			}
			return(.self);
		},
		get = function(i) {
			return(base::get(paste(i),dataEnv));
		},
		list = function() {
			if(n==0)	return(list());
			result = vector("list",n);
			for( i in 1:n) {
				result[[i]] = .self$get(i);
			}
			return(result);
		},
		unlist = function() {
			return(base::unlist(.self$list(), recursive=FALSE, use.names = FALSE));
		},
		show = function() {
			cat(".listBuilder1 object.\nIternal object in MatrixEQTL package.\n");
			cat("Number of elements:", .self$n, "\n");
		}
	)
)

.SlicedDataOffsets = function(x) {
	rez = integer(x$nSlices()+1);
	for( i in 1:x$nSlices()) {
		rez[i+1] = rez[i] + NROW( x$getSliceRaw(i) );
	}
	return(rez);
}


CorrMeta_testAll = function(
	snps1, gene1, cvrt1,
	snps2, gene2, cvrt2,
	pvThreshold) {

	gene.offsets = .SlicedDataOffsets(gene1);
	snps.offsets = .SlicedDataOffsets(snps1);
	
	gene.names = rownames(gene1);
	snps.names = rownames(snps1);
	
	### Test all the conditions
	{
		stopifnot(all( gene.names == rownames(gene2) ))
		stopifnot(all( snps.names == rownames(snps2) ))
		
		stopifnot(all( gene.offsets == .SlicedDataOffsets(gene2) ))
		stopifnot(all( snps.offsets == .SlicedDataOffsets(snps2) ))
		
		stopifnot( ncol(gene1) == ncol(snps1) );
		stopifnot( ncol(gene2) == ncol(snps2) );
		
		if(nrow(cvrt1) > 0)
			stopifnot( ncol(cvrt1) == ncol(snps1) );
			
		if(nrow(cvrt2) > 0)
			stopifnot( ncol(cvrt2) == ncol(snps2) );
	}
	
	ngene = tail(gene.offsets,1);
	nsnps = tail(snps.offsets,1);
	nsamples = ncol(gene1) + ncol(gene2);
	dfFull = nsamples - 2 - nrow(cvrt1) - nrow(cvrt2);
# 	crthreshfun = function(pv) {
# 		thr = qt(pv/2, dfFull, lower.tail = FALSE);
# 		thr = thr^2;
# 		thr = sqrt(  thr / (dfFull + thr) );
# 		thr[pv >= 1] = 0;
# 		thr[pv <= 0] = 1;
# 		return( thr );
# 	}
# 	crthreshold = crthreshfun(1e-5);
	ttthreshfun = function(pv) {
		thr = qt(pv/2, dfFull, lower.tail = FALSE);
		return( thr );
	}
	ttthreshold = ttthreshfun(pvThreshold);
	pvfun = function(x) { return( (pt(-abs(x),dfFull)*2)); }

	cat('Preprocessing data set 1','\n');
	worker1 = .SingleMatrixEQTL$new(snps1, gene1, cvrt1);
	cat('Preprocessing data set 2','\n');
	worker2 = .SingleMatrixEQTL$new(snps2, gene2, cvrt2);

	cat('Running analysis','\n');
	{
		collect.geneid = new('.listBuilder1');
		collect.snpsid = new('.listBuilder1');
		collect.tvalue = new('.listBuilder1');
	
		tic = proc.time();
	
		for( gsl in 1:gene1$nSlices() ) { # gsl = 1
			# cat( ssl, 'of', snps1$nSlices(), '\n');
		 	
			tt1 = worker1$get_tstats(gsl);
			tt2 = worker2$get_tstats(gsl);
			
			sum1 = sum2 = sum11 = sum22 = sum12 = 0;
			for( ssl in seq_along(tt1) ) { # ssl = 1
				sum1  = sum1  + rowSums(tt1[[ssl]])
				sum2  = sum2  + rowSums(tt2[[ssl]])
				sum11 = sum11 + rowSums(tt1[[ssl]]^2)
				sum22 = sum22 + rowSums(tt2[[ssl]]^2)
				sum12 = sum12 + rowSums(tt1[[ssl]]*tt2[[ssl]])
			}
			crosscor = (sum12 - sum1*sum2/nsnps) / sqrt( ( sum11 - sum1^2/nsnps)*( sum22 - sum2^2/nsnps) );
			
			### Rows - genes (slice gsl)
			### Columns - snps (slice ssl)

			for( ssl in seq_along(tt1)) { # sl=1
				tt = (tt1[[ssl]]+tt2[[ssl]])/sqrt(2+2*crosscor);
				selected = which(abs(tt) > ttthreshold, arr.ind = TRUE, useNames=FALSE);
				if(NROW(selected)>0) {
					collect.geneid$add(selected[,1] + gene.offsets[ssl]);
					collect.snpsid$add(selected[,2] + snps.offsets[ssl]);
					collect.tvalue$add(tt[selected]);
				}
			}
			
			cat(floor(100*gene.offsets[gsl+1]/tail(gene.offsets,1)),'% done, ',.s(collect.tvalue$count),' eQTLs found','\n',sep = '');
			flush.console()
			rm(tt1,tt2,tt,selected);
			# gc();
		}
	
		toc = proc.time();
		
		cat('Analysis finished in',(toc-tic)[3],'seconds','\n');
	}
	
	gene.factor = collect.geneid$unlist();
	levels(gene.factor) = gene.names;
	class(gene.factor) = 'factor'
	
	snps.factor = collect.snpsid$unlist();
	levels(snps.factor) = snps.names;
	class(snps.factor) = 'factor'
	
	result = data.frame( gene = gene.factor, snps = snps.factor, tstat = collect.tvalue$unlist(), 
								check.rows=FALSE, stringsAsFactors=FALSE);
	return(result);
}	
