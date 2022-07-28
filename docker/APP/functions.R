# order backbone points

order_base = function(backbone_points){
  backbone_points = backbone_points
  backbone_points$base_order = c(nrow(backbone_points):1)
  backbone_points = backbone_points[order(backbone_points$base_order,decreasing = FALSE),]
  
  f = vector()
  f[1] = 0
  
  sum = 0
  ## calculate distance between each backbone points
  for(i in 2:nrow(backbone_points)){
    sum = sum + sqrt((backbone_points$x[i]-backbone_points$x[i-1])^2 + (backbone_points$y[i]-backbone_points$y[i-1])^2)
    f[i] = sum 
  }
  
  
  backbone_points$length = f
  return(backbone_points)
}


# project points to base

project_points2base = function(backbone_points,query_points){
  base_points = backbone_points
  query_points = query_points[,c("x","y","z","distance2center","pos")]

  output = query_points[,c("x","y","z","distance2center","pos")]
  output$nn_index = 0
  output$nn_dist = 0
  output$nearest_Row = 0
  output$nearest_Column = 0
  output$shortest_path_order = 0
  output$shortest_path_length = 0
  output$note = "NA"
  
  for(i in c(1:nrow(query_points))){

    if(query_points[i,"pos"] %in% base_points$pos){
      output[i,"nn_index"] = -3
      output[i,"nn_dist"] = -3
      output[i,"nearest_Row"] = -3
      output[i,"nearest_Column"] = -3
      output[i,"shortest_path_order"] = -3
      output[i,"shortest_path_length"] = -3
      output[i,"note"] = "Base layer"
    }else{
      ### only project query points to base layer with larger radius
      f1 = base_points %>% subset(., distance2center >= query_points[i,"distance2center"])
      if(nrow(f1) == 0){
        output[i,"nn_index"] = -1
        output[i,"nn_dist"] = -1
        output[i,"nearest_Row"] = -1
        output[i,"nearest_Column"] = -1
        output[i,"shortest_path_order"] = -1
        output[i,"shortest_path_length"] = -1
        output[i,"note"] = "Query point with too large radius"
      }else{
        kd_closest = RANN::nn2(f1[,c('x','y')],query_points[i,c('x','y')],k=1,searchtype = "radius",radius = 5000)
        if(kd_closest$nn.idx == 0){
          output[i,"nn_index"] = -2
          output[i,"nn_dist"] = -2
          output[i,"nearest_Row"] = -2
          output[i,"nearest_Column"] = -2
          output[i,"shortest_path_order"] = -2
          output[i,"shortest_path_length"] = -2
          output[i,"note"] = "rann::nn2 failed"
        }else{
          output[i,"nn_index"] = kd_closest$nn.idx
          output[i,"nn_dist"] = kd_closest$nn.dists
          output[i,"nearest_Row"] = f1[kd_closest$nn.idx,c('x')]
          output[i,"nearest_Column"] = f1[kd_closest$nn.idx,c('y')]
          output[i,"shortest_path_order"] = f1[kd_closest$nn.idx,c('base_order')]
          output[i,"shortest_path_length"] = f1[kd_closest$nn.idx,c('length')]
          output[i,"note"] = "Successfully projected"
        }
      }
    }
    
    
  }
  return(output)
}

# calculate upper outlier per backbone point
zscore_per_backbone_point = function(converted_image,backbone_points){
  converted_image = converted_image

  backbone_points = backbone_points
  backbone_points$thickness_mean = 0
  backbone_points$thickness_sd = 0
  
  
  
  ## set outlier range per backbone base
  for(j in 1:nrow(backbone_points)){
    all_points_projected_here = subset(converted_image,nearest_Row == backbone_points[j,"x"] & nearest_Column == backbone_points[j,"y"])
    if(nrow(all_points_projected_here)==1){
      backbone_points[j,"thickness_mean"] = mean(all_points_projected_here$nn_dist)
      backbone_points[j,"thickness_sd"] = 0
    }else{
      backbone_points[j,"thickness_mean"] = mean(all_points_projected_here$nn_dist)
      backbone_points[j,"thickness_sd"] = sd(all_points_projected_here$nn_dist)
    }
    
  }
  
  return(backbone_points)
}



# filtering by angle
qc_angle_outlier = function(converted_image,x0,y0,backbone_points){
  x0 = x0
  y0 = y0
  
  converted_image = converted_image
  converted_image$angleCBS = 0
  converted_image$zscore = 0

  for(i in 1:nrow(converted_image)){
    
    
    if(converted_image[i,"note"] != "Successfully projected"){
      converted_image[i,"angleCBS"] = "Projection failed"
      converted_image[i,"zscore"] = "Projection failed"
    }else{
      ## angle
      A=c(x0,y0)
      B=c(converted_image$nearest_Row[i],converted_image$nearest_Column[i])
      C=c(converted_image$x[i],converted_image$y[i])
      vector1 = c(A[1] - B[1], A[2] - B[2])
      vector2 = c(C[1] - B[1], C[2] - B[2])
      num = (vector1[1] * vector2[1] + vector1[2] * vector2[2])
      den = sqrt(vector1[1]^2 + vector1[2]^2) * sqrt(vector2[1]^2 + vector2[2]^2)
      angle = acos(num/den)
      angle = (360 * angle)/(2 * pi)
      converted_image[i,"angleCBS"] = angle
      
      ## outlier
      outlier_value = subset(backbone_points,x == converted_image[i,"nearest_Row"] & y == converted_image[i,"nearest_Column"])
      # print(outlier_value)
      
      if(outlier_value[1,"thickness_sd"] == 0){
        converted_image[i,"zscore"] = 0
      }else{
        converted_image[i,"zscore"] = (converted_image[i,"nn_dist"]-outlier_value[1,"thickness_mean"])/outlier_value[1,"thickness_sd"]
      }
      
      
    }
    
  }

  return(converted_image)
  
}

# scale marker to 0-1
scale_marker = function(x){(x-min(x))/(max(x)-min(x))}
