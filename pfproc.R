library(hash)
library(reticulate)
library(installr)
library(splines)
library(httr)
library(future)
library(jsonlite)
library(dequer)
library(R6)

#flsa<-import_from_path('flsa', path="C:\\Users\\zulfi\\OneDrive\\python\\convex")


evaluate_y_axis<-function( q, x ){
  y<-interp1(q.keys(),q.values(), x, type="cubic")
  y
}

# Return a spline from xy hash

smoothen_scaled_histogram<-function( xy_hash ){
    kx<-keys(xy_hash)
    x<-as.numeric(kx)
    if (length(x)<1){
      return(NULL)
    }
    #print(paste('in smoothen_scaled_histogram:', length(kx)))
    mx<-max(x)
    mn<-min(x)
    y<-c()
    for (i in kx){
      tryCatch(
        {
          #print(xy_hash[[i]])
          v<-as.numeric(xy_hash[[i]])
          #print(paste('k=',i,'; v=',v))
          y<-append(y,v)          
        }, 
        error=function(e){e}
      )
    }
    #print(paste('len x:',length(x)))
    xp<-(as.vector(x)-mn)/(mx-mn)
    ans<-smooth.spline( xp, y, df=16)
    z<-predict(ans,xp)$y
    z<-z/sum(as.vector(z))
    ans<-smooth.spline( xp, z, df=16)
    ans
}


kullback_leibler_divergence_hash<-function(q1,q2){
  xaxis<-c(q1$x,q2$x)
  xp<-xaxis[1]
  ans<-0
  for (xv in xaxis){
    q1v<-evaluate_y_axis( q1, xv)
    q2v<-evaluate_y_axis( q2, xv)
    dx<-xv-xp
    dans<-dx*( q1v*log( q1v/(q2v+0.00001)))
    ans<-ans+dans
    xp<-xv
  }
  ans
}

kullback_leibler_divergence<-function(q1,q2){
  x<-seq(0,1,by=0.01)
  ans<-0
  tryCatch({
    if ( !is.null(q1)  && !(is.null(q2))) {
      y1<-as.vector(predict( q1, x)$y)
      y2<-as.vector(predict( q2, x)$y)
      ans<-sum( 0.01*y1*log(((y1+0.000001)/(y2+0.0001)+0.0001)))
    }
  }, error=function(e){e})
  return(ans)
}
  
traverse_count<-function( node, H ){
  ht<-as.character(attr(node, "height"))
  if (length(ht)>0){
    if (has.key(ht,H)){
      v<-H[[ht]]
      tryCatch({
        H[[ht]]<-as.numeric(v)+1
      },error=function(e){e})
    } else {
      H[[ht]]<-1
    }
  }
  n<-as.numeric(attr(node, "members"))
  if ((!is.empty(n)) && (n>0) ){
    for (childnode in node){
      traverse_count( childnode, H)
    }
  }
  H
}


tree_probability_density<-function( tree ){
  # traverse the tree filling up a key-value dictionary
  # where value is just the number of members
  count_by_height<-hash()
  clear(count_by_height)
  count_by_height<-traverse_count( tree, count_by_height )
  probability_density<-smoothen_scaled_histogram( count_by_height )
  probability_density
}

check_tree<-function( tree ){
  q<-tree_probability_density(tree)
  if (is.null(q)){
    print('probability density is empty')
    return(FALSE)
  }
  #print('Is spline q working ok:')
  ans<-predict( q, c(0.2,0.4, 0.6, 0.8))

  #print(ans$y)
  if (length(ans$y)==4){
    #print("everything is fine")
    return(TRUE)
  } else {
    print("nope")
    return(FALSE)
  }
}
dendrogram_distance<-function( T1, T2 )
{
  q1<-tree_probability_density( T1 )
  q2<-tree_probability_density( T2 )
  
  #these should be spline objects
  #print('Is spline q1 working ok:')
  #print(predict( q1, c(0.2,0.4, 0.6, 0.8)))
  #print('Is spline q2 working ok:')
  #print(predict( q2, c(0.2,0.4, 0.6,0.8)))
  
  qdist<-abs(kullback_leibler_divergence( q1, q2 ))
  if (is.na(qdist)){
    qdist<-1e-4
  }
  qdist
}

split_dataset<-function( D, size=500 ){
  N<-dim(D)[1]
  chunks<-list()
  ptr<-1
  while (ptr<N){
    chunk<-D[ptr:(ptr+size-1),]
    chunk[is.na(chunk)]<-0
    chunks<-append(chunks, list(chunk))
    ptr<-ptr+size
  }
  chunks
}

hierarchical_cluster<-function(data){
  dissimilarity<-dist( data, method="euclidean")
  H<-hclust( dissimilarity, method="ward.D2")
  H
}

async_hierarchical_cluster<-function(data){
  # let's assume blocking response for now
  print(paste('first line', data[1,]))
  res<-POST( url='http://127.0.0.1:6062/', body=toJSON(data))
  #headers<-headers(res)
  #print(headers)
  res_json<-content(res,as="text")
  
  x<-fromJSON(res_json)
  tree<-unserializeJSON(res_json)
  tree
}

Unit_test_serialize_unserialize_json<-function(){

  random_data<-data.frame( x=rnorm(30), y=rnorm(30), z=rnorm(30))
  dissimilarity<-dist( random_data, method="euclidean")
  tree1<-hclust( dissimilarity )
  tree2<-unserializeJSON(serializeJSON(tree1))
  all.equal( tree1, tree2)
}

convex_cluster<-function( data ){
   convex_tree_data<-flsa$flsa( data )
   convex_tree_data
}


serial_cluster_medium_data<-function( mediumDataset ){
  epsilon<-5e-6
  chunks<-split_dataset( mediumDataset )
  trees<-c()
  bigtree<-NULL
  prev_bigtree<-NULL
  track_small_increment<-0
  num_small_inc_for_stop<-5
  iteration<-0
  for ( chunk in chunks ){
    tree<-hierarchical_cluster( chunk )
    if (is.null(bigtree)){
      bigtree<-tree
    } else {
      prev_bigtree<-bigtree
      bigtree<-merge(as.dendrogram(bigtree), as.dendrogram(tree))
    }
    if (! is.null( prev_bigtree) ){
      if ( check_tree( prev_bigtree) &&
           check_tree( bigtree) ){
        
        cur_dist<-dendrogram_distance( bigtree, prev_bigtree)
      
        print(paste("iteration:", iteration,
                    "tree dist increment:", cur_dist))
        if (cur_dist < epsilon){
          track_small_increment<-track_small_increment+1
        } else {
          track_small_increment<-0
        }
        if ( track_small_increment > num_small_inc_for_stop){
        # stop the process
          print( 'stopping the clustering')
        }
        iteration<-iteration+1
      }
    }
  }
}


parallel_cluster_medium_data<-function( mediumDataset ){
  plan(sequential)
  epsilon<-5e-6
  chunks<-split_dataset( mediumDataset )
  n_chunks<-length(chunks)
  
  bigtree<-NULL
  prev_bigtree<-NULL
  track_small_increment<-0
  num_small_inc_for_stop<-5
  iteration<-0
  # make this part asynchronous
  task_list<-c()
  tree_list<-vector("list", n_chunks)

  for ( job_id in 1:n_chunks ){
    f<-future({async_hierarchical_cluster( chunks[[job_id]] )})
    task_list[[job_id]]<-f
  }
  
  n_tasks<-n_chunks
  print(n_tasks)
  
  # wait for any jobs
  
  jobs_done<-rep(0,n_tasks)
  still_working<-TRUE
  
  while (still_working){

    for (jj in 1:n_tasks){
      if ((jobs_done[jj]==0) && resolved(task_list[[jj]])) {
        print( paste('completed job number ', jj, ' out of ', n_tasks))
        
        tree<-as.dendrogram(value(task_list[[jj]]))
        tree_list[[jj]]<-tree
        
        jobs_done[jj]<-1
        
        print(paste('extracted tree ', jj, ' of class ', class(tree)))
        
        if (is.null(bigtree)){
          bigtree<-tree
        } else {
          prev_bigtree<-bigtree
          tryCatch({
            print( attr(as.dendrogram(bigtree),"members") )
            bigtree<-merge(prev_bigtree,tree)
          
          }, error=function(e){e})
          #bigtree<-merge(as.dendrogram(bigtree), as.dendrogram(tree))
        }
        
        print(paste('bigtree size is now', attr(as.dendrogram(bigtree), "members")))
        
        if (! is.null( prev_bigtree) ){
          
          #check_pt <- check_tree( prev_bigtree)
          #check_t <- check_tree( bigtree)
          
          check_pt<-TRUE
          check_t<-TRUE
          
          if (! check_pt ){
            print( 'validation failed on prev big tree')
          }
          
          if (! check_t ){
            print( 'validation failed on big tree')
          }
          
          if ( check_pt  && check_t ){
        
            cur_dist<-dendrogram_distance( bigtree, prev_bigtree)
        
            print(paste("iteration:", iteration,
                      "tree dist increment:", cur_dist))
            if (cur_dist < epsilon){
              track_small_increment<-track_small_increment+1
            } else {
              track_small_increment<-0
            }
            if ( track_small_increment > num_small_inc_for_stop){
            # stop the process
              print( 'stopping the clustering')
            }
            iteration<-iteration+1
          }
        }
      }
    }
    print(jobs_done)
    print(all(jobs_done==1))
    if (all(jobs_done==1)){
      still_working<-FALSE
    }
  }
  as.dendrogram(bigtree)
}


proc<-R6Class( "Process", list(

    epsilon = 5e-6,
    current_error = 1000,
    chunks = NULL,
    prev_bigtree = NULL,
    bigtree = NULL, 
    task_list = NULL,
    tree_list = NULL,
    jobs_processed = NULL,
    jobs_dispatched = NULL,
    jobs_queue = NULL,
    errors = NULL,

    initialize = function( dataset ) {
      self$chunks<-split_dataset( dataset )
      n<-length(self$chunks)
      self$task_list<-vector( "list", n )
      self$tree_list<-vector( "list", n )
      self$jobs_processed<-rep(0,n)
      self$jobs_dispatched<-rep(0,n)
      self$errors<-c()

      self$load_jobs( 5 )
    },

    load_jobs= function(k){
      m<-1
      for ( rr in 1:length(self$chunks)){
        if ( !resolved(self$task_list[[rr]])) {
          if ( m <= k){
            pushback( self$jobs_queue, rr)
            m<-m+1
          }
        } 
      }
    },

    dispatch_jobs = function( ){
      while ( length(self$jobs_queue) > 0 ){
        id = pop( self$jobs_queue)
        f<-future({async_hierarchical_cluster( self$chunks[[id]] )})
        self$task_list[[id]]<-f
        self$jobs_dispatched[id]<-1
      }    
   },

    print_message = function( ){
      k<-length(errors)
      print( paste( k, errors[k] ))
    }

    handle_job_done = function( id ){
      
      if ( ! resolved( self$task_list[[id]] ) ){
        return(FALSE)
      }
      self$jobs_processed[id]<-1

      self$tree_list[[id]] <- as.dendogram( value( self$task_list[[ id ]] ) )
      tree <- self$tree_list[[ id ]]

      self$prev_bigtree <- self$bigtree
      self$bigtree <- merge( prev_bigtree, self$bigtree )

      self$current_error<-self$dendogram_distance( self$bigtree, self$prev_bigtree)
      self$errors<-append( self$errors, self$current_error )
      self$print_message()      
    },
    
    check_termination_condition = function( ) {
      if (sum(self$jobs_processed) > length(self$chunks)-1 ){
        return(TRUE)
      }
      if ( self$current_error > 0 && self$cureent_error < epsilon ){
        return(TRUE)
      }
      return(FALSE)
    },

    run = function(){
      while (TRUE){
        finish<-self$check_termination_condition()
        if ( finish ) {
          return( self$bigtree )
        }

        # send out at most 5 more jobs if processed jobs
        # exceed done jobs - 5
        if ( sum(self$jobs_processed) + 5 >= sum(self$jobs_dispatched) ){
          # select and dispatch at most five more jobs
          self$load_jobs(5)
          self$dispatch_jobs()
        }

        # look through jobs done
        for ( id in 1:length(self$chunks)){
          self$handle_job_done(id)
        }
      }
      return self$bigtree
    }
  )
)


parallel_controlled_cluster_medium_data<-function( mediumDataset ){
  P<-Process()
  P$initialize( mediumDataset)
  bt<-P$run()
}

# Old code
parallel_controlled_cluster_medium_data_old<-function( mediumDataset ){
  plan(sequential)
  epsilon<-5e-6
  current_error<-1000
  chunks<-split_dataset( mediumDataset )
  n_chunks<-length(chunks)
  
  bigtree<-NULL
  prev_bigtree<-NULL
  track_small_increment<-0
  num_small_inc_for_stop<-5
  iteration<-0

  # make this part asynchronous
  
  jobs_queue<-queue()
  last_job<-1
  
  load_jobs<-function( ){
    for (k in 1:5 ){
      pushback( jobs_queue, last_job+k)
    }
  }
  
  send_job<-function( id ){
    f<-future({async_hierarchical_cluster( chunks[[id]] )})
    task_list[[id]]<-f    
  }
  
  send_jobs_wave<-function(){
    while (length(jobs_queue)){
      send_job( pop(jobs_queue))
    }
  }
  
  task_list<-vector("list", n_chunks)
  tree_list<-vector("list", n_chunks)

  # wait for any jobs
  
  jobs_done<-rep(0,n_chunks)
  still_working<-TRUE
  
  while (still_working){
    
    if (last_job< 5){
      load_jobs()
      send_jobs_wave()
    }
    
    if (last_job> 5){
      if (sum(jobs_done)+5 > last_job){
        load_jobs()
        send_jobs_wave()
      }
    }
    
    for (jj in 1:n_chunks){
      if ((jobs_done[jj]==0) && (!is.null(task_list[[jj]])) && resolved(task_list[[jj]])) {
        print( paste('completed job number ', jj, ' out of ', n_chunks))
        
        tree<-as.dendrogram(value(task_list[[jj]]))
        tree_list[[jj]]<-tree
        
        jobs_done[jj]<-1
        
        print(paste('extracted tree ', jj, ' of class ', class(tree)))
        
        if (is.null(bigtree)){
          bigtree<-tree
        } else {
          prev_bigtree<-bigtree
          tryCatch({
            print( attr(as.dendrogram(bigtree),"members") )
            bigtree<-merge(prev_bigtree,tree)
            
          }, error=function(e){e})
          #bigtree<-merge(as.dendrogram(bigtree), as.dendrogram(tree))
        }
        
        print(paste('bigtree size is now', attr(as.dendrogram(bigtree), "members")))
        
        if (! is.null( prev_bigtree) ){
          
          #check_pt <- check_tree( prev_bigtree)
          #check_t <- check_tree( bigtree)
          
          check_pt<-TRUE
          check_t<-TRUE
          
          if (! check_pt ){
            print( 'validation failed on prev big tree')
          }
          
          if (! check_t ){
            print( 'validation failed on big tree')
          }
          
          if ( check_pt  && check_t ){
            
            cur_dist<-dendrogram_distance( bigtree, prev_bigtree)
            
            print(paste("iteration:", iteration,
                        "tree dist increment:", cur_dist))
            if (cur_dist>0){
              current_error<-cur_dist
            }
            if (cur_dist < epsilon){
              track_small_increment<-track_small_increment+1
            } else {
              track_small_increment<-0
            }
            if ( track_small_increment > num_small_inc_for_stop){
              # stop the process
              print( 'stopping the clustering')
            }
            iteration<-iteration+1
          }
        }
      }
    }

    if ( current_error< epsilon){
      still_working<-FALSE
    }
    if (all(jobs_done==1)){
      still_working<-FALSE
    }
  }
  as.dendrogram(bigtree)
}














######################

cluster_very_large_data<-function( largeDataset ){
  m<-largeDataset/500
  epsilon<-0.01
  error<-100.
  prev_tree<-c()
  current_tree<-c()
  while (error> epsilon){
    if (jobs_done >= m){
      break
    }
    next_job<-next_job_for_submission( )
    data_chunk<-largeDataset[ (500*next_job):(501*next_job),]
    submit( next_job, data_chunk)
    tree_array<check_jobs_results()
    if (tree_array$n>=1){
      prev_tree<-current_tree
      for (t in tree_array$trees){
        current_tree<-merge(current_tree, t)
      }
      error<-dendrogram_distance( current_tree, prev_tree)
      print(paste('error:', error))
    }
    if (error<epsilon){
      break
    }
  }
}
