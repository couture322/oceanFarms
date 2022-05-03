### Messing with K


sim_fisheryAqK <-
  function(fish,
           fleet,
           manager,
           num_patches = 10,
           sim_years = 1,
           burn_years = 25,
           crashed_pop = 1e-3,
           random_mpas = F,
           farmSize=0, # to set regular MPA/farms of a specific size (value is numer of patches, will be distributed evenly along the total area)
           enviro = NA,
           enviro_strength = 1,
           rec_driver = "stochastic",
           est_msy = F,
           time_step,
           max_window = 10,
           min_size = 1,
           mpa_habfactor = 1,
           farm_yrs = NA, # keep as NA for normal mpa, add vector for variable farm presence after mpa year. Length should be length(mpaYr:simYr)
           fallowFactor = 1,
           farmAttr=1,#value of 1 does nothing, higher increases movement to farm
           ddAttr=1,
           #bufferMove=1,
           buffSize=2,
           mImp=1, #impact of the farm on natural mortality rate: <1 improves health, >1 damaging
           sprinkler = FALSE,
           keep_burn = FALSE,
           tune_costs = FALSE,
           patchK = 10000000000,  # carrying capacity per patch (in biomass)
           prodIncr=1, #productivity increase with retention
           lag = 0) { # when does retention start? right away? (lag = 0) or after x # of years
    msy <- NA
    
    p_msy <- NA
    
    e_msy <- NA
    
    max_r_msy <-  NA
    
    b0 <- NA
    
    farmYrs<-rep(FALSE,sim_years)
    
    if(is.na(farm_yrs))
      {
        farmYrs <- rep(TRUE,length(farmYrs))
        } else {
          farmYrs[manager$year_mpa:sim_years] <-farm_yrs
    }
    
    ddwaa <- data.frame(waa=c(0,fish$weight_at_age))%>%
      mutate(age=1:length(waa),
             statGrwth=c(waa[1],waa[2:length(waa)]-lag(waa)[2:length(waa)]))
    
    if (sprinkler == FALSE & mpa_habfactor == 1){
      burn_years <- 1
    }
    
    if (est_msy == T) {
      e_msy_guess <- fish$m / fleet$q
      
      tol <- 100
      
      lower <- 0
      
      upper <- e_msy_guess * 3
      
      golden <- (sqrt(5) - 1) / 2
      
      best <- 1000
      
      delta_best <- 100
      
      counter <-  0
      
      while (delta_best > tol | counter < 20) {
        counter <- counter + 1
        
        constant <- (1 - golden) * (upper - lower)
        
        x1 <- lower + constant
        
        x2 <- upper - constant
        
        yield_1 <-
          estimate_msy(x1,
                       fish = fish,
                       fleet = fleet,
                       num_patches = 1)
        
        yield_2 <-
          estimate_msy(x2,
                       fish = fish,
                       fleet = fleet,
                       num_patches = 1)
        
        delta_best <-  (best -  min(yield_1, yield_2)) ^ 2
        
        best <- min(yield_1, yield_2)
        
        if (yield_1 < yield_2) {
          lower <- lower
          
          upper <- x2
        } else{
          lower <- x1
          
          upper <- upper
          
        }
        
        
      } # close golden while
      
      
      msy_foo <- function(seed, lower, upper, fish, fleet) {
        msy_fit <-
          nlminb(
            mean(c(lower, upper)),
            estimate_msy,
            fish = fish,
            fleet = fleet,
            lower = 0,
            seed = seed,
            num_patches = 1
          )
        
        out <- list()
        
        out$e_msy <- msy_fit$par
        
        out$msy <- -msy_fit$objective
        
        ref_points <-
          estimate_msy(
            out$e_msy,
            fish = fish,
            fleet = fleet,
            use = "other",
            seed = seed
          )
        
        out$b_msy <- ref_points$b_msy
        
        out$r_msy <- ref_points$r_msy
        
        return(out)
        
      } # close msy foo
      
      msy_ests <-
        map(
          sample(1:10000, 1, replace = F),
          msy_foo,
          fish = fish,
          fleet = fleet,
          lower = lower,
          upper = upper
        )
      
      max_r_msy <- max(map_dbl(msy_ests, "r_msy"))
      
      e_msy <- mean(map_dbl(msy_ests, "e_msy"))
      
      fleet$e_msy <- e_msy
      
      msy <- mean(map_dbl(msy_ests, "msy"))
      
      fish$msy <- msy
      
    } # close if estimate msy
    
    
    sim_years <- burn_years + sim_years
    
    if (fleet$fleet_model == "supplied-catch") {
      sim_years <- sim_years + 1
    }
    
    
    pop <-
      expand.grid(
        year = 1:sim_years,
        patch = 1:num_patches,
        age = seq(fish$min_age, fish$max_age, fish$time_step)
      ) %>%
      dplyr::mutate(
        numbers = 0,
        biomass = 0,
        ssb = 0,
        numbers_caught = 0,
        biomass_caught=0,
        profits = 0,
        effort = 0,
        f = 0,
        mpa = F,
        cost = 0
      ) %>%
      dplyr::as_data_frame() %>%
      arrange(year, patch, age)
    
    
    effort <- vector(mode = "double", length = sim_years)
    
    f <- vector(mode = "double", length = sim_years)
    
    
    if (rec_driver == "stochastic") {
      # rec_devs <-
      #   rnorm(
      #     sim_years,
      #     mean = -(fish$sigma_r ^ 2) / 2,
      #     sd = fish$sigma_r
      #   )
      
      rec_devs <-
        rnorm(sim_years,
              mean = 0,
              sd = fish$sigma_r)
      
      ## autocorrelated recruitment deviations
      for (t in 2:length(rec_devs)) {
        rec_devs[t] <-
          rec_devs[t - 1] * fish$rec_ac + sqrt(1 - fish$rec_ac ^ 2) * rec_devs[t]
      }
    } else if (rec_driver == "environment") {
      if (length(enviro) != sim_years) {
        stop("environment must be same length as sim_years")
      }
      
      rec_devs <-
        rnorm(sim_years,
              mean = enviro_strength * enviro,
              sd = fish$sigma_r)
    }
    
    effort_devs <-
      rnorm(sim_years,
            mean = 0,
            sd = fleet$sigma_effort)
    
    for (t in 2:length(effort_devs)) {
      effort_devs[t] <-
        effort_devs[t - 1] * fleet$effort_ac + sqrt(1 - fleet$effort_ac ^ 2) * effort_devs[t]
    }
    
    # mpa_locations <- -1
    
    prop_mpas <- round(num_patches * manager$mpa_size)

  if(prop_mpas==farmSize){
      mpa_locations <- c(1:farmSize)
    
      } else{
      
        if(farmSize>0 & manager$mpa_size>0) {
      
        nFrms<-round(prop_mpas/farmSize)
      
        farmStart<-round(seq(0,num_patches,length.out = nFrms+2))[0:nFrms]+1
      
        farmLFunc<-function(x){
        seq(x,(x+farmSize-1),by=1)
      
          }
      
      locsList<-lapply(farmStart, farmLFunc)
      mpa_locations<-unlist(locsList)
      
    } else {
      
      if (random_mpas == T & prop_mpas > 0) {
        
        ms <- min(prop_mpas, max(1, min_size * num_patches)) # calculate min MPA size for random MPA generation
        
        cwidth <- num_patches / ms
        
        atemp <- tibble(patch = 1:num_patches) %>%
          dplyr::mutate(cluster = cut(patch, pmax(2,round(cwidth))))
        btemp <-
          sampling::cluster(atemp,
                            cluster = "cluster",
                            pmin(n_distinct(atemp$cluster),ceiling(prop_mpas / ms)),
                            method = "srswor")
        
        ctemp <- sampling::getdata(atemp, btemp) %>%
          sample_n(pmin(prop_mpas, nrow(.)))
        
        mpa_locations <- ctemp$patch
        
        
      } else {
        mpa_locations <-
          (1:num_patches)[0:prop_mpas] #weird zero is in case prop_mpas is zero
      }
      
      if (!all(is.na(manager$mpa_locations))){
        
        if (prop_mpas > 0){
          warning("overwriting MPA size with specific MPA locations")
        }
        
        mpa_locations <- manager$mpa_locations
        
        if (max(mpa_locations) > num_patches){
          stop("invalid MPA location supplied, make sure MPAs fit inside number of patches")
        }
      }
    }
  }
    
    habitat <- rep(1, num_patches)
    habitat[mpa_locations]  <- mpa_habfactor # define habitat factor in MPA other than 1 (better (>1) or worse? (<1))
    
    fallowHab <- rep(1,num_patches)
    fallowHab[mpa_locations] <- fallowFactor
    
    if (prop_mpas == 0){ # in case you want to simulate MPA habitat but without the MPAs
      
      mpa_locations <- (1:num_patches)[0:prop_mpas] #weird zero is in case prop_mpas is zero
      
    }
    
    ### buffer code
    
    bufferLocs<-if(exists("farmStart",envir = globalenv())==T){
      
      c(farmStart-buffSize,farmStart+farmSize+buffSize)
      
    }
    
    n0_at_age <-
      (fish$r0 / num_patches) * exp(-fish$m * seq(fish$min_age, fish$max_age, fish$time_step))
    
    n0_at_age[fish$max_age + 1] <-
      n0_at_age[fish$max_age + 1] / (1 - exp(-fish$m))
    
    b0_at_age <- n0_at_age * fish$weight_at_age
    
    ssb0_at_age <- n0_at_age * fish$ssb_at_age
    
    # generate time series of price, cost, and q if called for
    
    price_series <-
      generate_timeseries(
        fish$price,
        cv = fish$price_cv,
        ac = fish$price_ac,
        percent_slope = fish$price_slope,
        time = sim_years
      )
    
    q <-
      generate_timeseries(
        fleet$q,
        cv = fleet$q_cv,
        ac = fleet$q_ac,
        percent_slope = fleet$q_slope,
        time = sim_years
      )
    
    if (length(q) == 1) {
      q <- rep(q, sim_years)
    }
    if (length(price_series) == 1) {
      price_series <- rep(price_series, sim_years)
    }
    
    # tune costs based on some heavy fishing at b0
    
    hyp_f <- fish$m #hypothetical f
    
    hyp_effort <- hyp_f / mean(q[(burn_years + 1):sim_years])
    
    hyp_f_at_age <- hyp_f * fleet$sel_at_age
    
    hyp_b0_catch <-
      sum((hyp_f_at_age / (hyp_f_at_age + fish$m))  * b0_at_age * (1 - exp(-(
        hyp_f_at_age + fish$m
      ))))
    
    b0_revenue <-
      mean(price_series[(burn_years + 1):sim_years]) * hyp_b0_catch
    
    hyp_profits_guess <- b0_revenue * (1 - fleet$max_cr_ratio)
    
    cost_guess <-
      (b0_revenue - hyp_profits_guess) / hyp_effort ^ fleet$beta
    
    fleet$theta <-
      (fleet$max_perc_change_f * hyp_effort) / (hyp_profits_guess / hyp_effort)
    
    cost_series <-
      generate_timeseries(
        cost_guess,
        cv = fleet$cost_cv,
        ac = fleet$cost_ac,
        percent_slope = fleet$cost_slope,
        time = sim_years
      )
    
    cost_series <- (cost_series / max(cost_series)) * cost_guess
    
    fleet$cost <- cost_guess
    
    price_frame <-
      dplyr::data_frame(year = 1:sim_years, price = price_series)
    
    cost_frame <- dplyr::data_frame(year = 1:sim_years, cost = cost_series)
    
    pop <- pop %>%
      dplyr::select(-cost) %>%
      dplyr::left_join(cost_frame, by = "year") %>%
      dplyr::left_join(price_frame, by = "year")
    
    pop$numbers[pop$year == 1] <- rep(n0_at_age, num_patches)
    
    if (fleet$cost_function == "distance from port") {
      cost_frame <-
        expand.grid(year = 1:sim_years, patch = 1:num_patches) %>%
        dplyr::as_data_frame() %>%
        dplyr::left_join(cost_frame, by = "year") %>%
        dplyr::mutate(cost = cost * (1 + fleet$cost_slope * (patch - 1)))
      
      pop <- pop %>%
        dplyr::select(-cost) %>%
        dplyr::left_join(cost_frame, by = c("patch", "year"))
    }
    
    
    pop <- pop %>%
      dplyr::left_join(
        dplyr::data_frame(
          age = seq(fish$min_age, fish$max_age, fish$time_step),
          ssb_at_age = fish$ssb_at_age,
          weight_at_age = fish$weight_at_age
        ),
        by = "age"
      )
    
    y <- 1
    
    model_phase <- "burn"
    
    # adult_move_grid <-
    #   expand.grid(
    #     source = 1:num_patches,
    #     sink = 1:num_patches
    #   ) %>%
    #   mutate(
    #     distance = source - sink,
    #     prob = 1 / ((2 * pi) ^ (1 / 2) * fish$adult_movement) * exp(-(distance) ^
    #       2 / (2 * fish$adult_movement ^ 2))
    #   ) %>%
    #   group_by(source) %>%
    #   mutate(prob_move = prob / sum(prob))
    
    adult_move_grid <-
      expand.grid(from = 1:num_patches, to = 1:num_patches) %>%
      as.data.frame() %>%
      dplyr::mutate(distance = purrr::map2_dbl(from, to, ~ min(
        c(abs(.x - .y),
          .x + num_patches - .y,
          num_patches - .x + .y)
      ))) %>%
      dplyr::mutate(movement = ifelse(
        is.finite(dnorm(distance, 0, fish$adult_movement)),
        dnorm(distance, 0, fish$adult_movement),
        1
      ))  %>%
      group_by(from) %>%
      dplyr::mutate(prob_move = movement / sum(movement))
    
    
    adult_move_matrix <- adult_move_grid %>%
      ungroup() %>%
      dplyr::select(from, to, prob_move) %>%
      spread(to, prob_move) %>%
      dplyr:: select(-from) %>%
      as.matrix()
    
    larval_move_grid <-
      expand.grid(from = 1:num_patches, to = 1:num_patches) %>%
      as.data.frame() %>%
      dplyr::mutate(distance = purrr::map2_dbl(from, to, ~ min(
        c(abs(.x - .y),
          .x + num_patches - .y,
          num_patches - .x + .y)
      ))) %>%
      dplyr::mutate(
        larval_movement = ifelse(
          from %in% mpa_locations &
            sprinkler != FALSE,
          fish$larval_movement * sprinkler,
          fish$larval_movement
        )
      )
    
    
    larval_move_grid <-  larval_move_grid %>%
      dplyr::mutate(movement = ifelse(is.finite(
        dnorm(distance, 0, fish$larval_movement)
      ), dnorm(distance, 0, larval_movement), 1))  %>%
      group_by(from) %>%
      dplyr::mutate(prob_move = movement / sum(movement))
    
    larval_move_matrix <- larval_move_grid %>%
      ungroup() %>%
      dplyr::select(from, to, prob_move) %>%
      spread(to, prob_move) %>%
      dplyr::select(-from) %>%
      as.matrix()
    
    
    ##### simulation (post-burn) #####

    farm_attr <- expand.grid(to=mpa_locations,distance=0:buffSize)%>% # limit attraction to areas close (within buffsize) to the farm (each patch of the farm)
      as.data.frame()%>%
      mutate(farmImpacts=farmAttr) # add the farmAttr value, only these will merge into the adult_move_grid df later with the appropriate distances
    
    
    for (y in 1:(sim_years - 1)) {
      # Move adults
      
      now_year <- pop$year == y
      
      if (y > burn_years) { #fish$density_movement_modifier < 1 & 
        
        slope <- fish$adult_movement - (fish$adult_movement * fish$density_movement_modifier)
        
        # depletion = seq(0,2, by = 0.1)
        #
        # pmin(fish$adult_movement, slope * depletion + (fish$adult_movement * fish$density_movement_modifier)) %>%
        #   plot()
        
        ### LAG
        
        # if(lag==0) {
        #  # notLeave <- get(ifelse(length(mpa_locations)>0,"mpaLv","noMpaLv"))
        #   farm_attr2<-farm_attr
        #   } else { 
        #     if(lag>0 & y < (manager$year_mpa + lag)){
        #     
        #      # notLeave<-noMpaLv
        #       farm_attr2<-farm_attr%>%mutate(farmImpacts=1)
        #     
        #     } else {
        #     #  notLeave<-mpaLv
        #       farm_attr2<-farm_attr
        #     }
        #   
        # }
        

        how_crowded <- pop %>%
          filter(now_year) %>%
          group_by(patch) %>%
          summarise(ssb = sum(ssb),.groups="drop") %>%
          arrange(patch) %>%
          mutate(depletion = ssb / fish$ssb0) %>%
          mutate(move_rate = pmin(
            fish$adult_movement,
            slope * depletion + (fish$adult_movement * fish$density_movement_modifier)
          )) %>%
          select(patch, depletion, move_rate)

        ##Denstity - dependent attraction
        # ddAttrDf<-how_crowded%>%
        #   select(-move_rate)%>%
        #   mutate(ddAtt = case_when(ddAttr == 1 ~ 0,
        #                            depletion >= 0.5 ~ (depletion * ddAttr),
        #                            TRUE ~ 0))%>%  # use depletion instead of pK since I don't want to change this relationship manually (by changing pK)
        #   select(to=patch,ddAtt)
        
        if(y > manager$year_mpa){
          
          ddAttrDf<-how_crowded%>%
            arrange(desc(depletion))%>%
            slice(1:(num_patches*manager$mpa_size))%>%
            filter(depletion > 0.3)%>%
            mutate(ddAtt=ddAttr)%>%
            select(to=patch,ddAtt,depletion)
        } else {
          
          ddAttrDf<-how_crowded%>%
            arrange(desc(depletion))%>%
            slice(1:(num_patches*manager$mpa_size))%>%
            mutate(ddAtt=0)%>%
            select(to=patch,ddAtt,depletion)
          
        }
          
        #print(ddAttrDf)
        
        adult_move_grid <-
          expand.grid(from = 1:num_patches, to = 1:num_patches) %>%
          as.data.frame() %>%
          left_join(how_crowded%>%
                      select(patch,move_rate), 
                    by = c("from" = "patch")) %>%
          dplyr::mutate(distance = purrr::map2_dbl(from, to, ~ min(
            c(abs(.x - .y),
              .x + num_patches - .y,
              num_patches - .x + .y)
          ))) %>%
          dplyr::mutate(movement = ifelse(
            is.finite(dnorm(distance, 0, move_rate)),
            dnorm(distance, 0, move_rate),
            1
          ))
        
        
        if((y-burn_years) >= manager$year_mpa){ # apply attractor values just after farm is deployed
          adult_move_grid <- adult_move_grid %>%
            ungroup()%>%
            left_join(.,farm_attr,by=c("to","distance"))%>% # merge attractor values with the 'to' patches
            left_join(.,ddAttrDf,by= "to")%>% # merge attractor values with the 'to' patches
            mutate(combAttr=replace_na(farmImpacts,1) + replace_na(ddAtt,0),
                   movement=movement*combAttr)%>%## adjust 'prob_move' by 'farm_attr'; increase probabilty that move FROM farm locations
            select(-combAttr)%>%
            group_by(from) %>%
            dplyr::mutate(prob_move = movement / sum(movement))#,
          #  prob_move = prob_move*noLve)
          
          
        } else {
          adult_move_grid <- adult_move_grid%>%
            group_by(from) %>%
            dplyr::mutate(prob_move = movement / sum(movement))
        }
        
        
        adult_move_matrix <- adult_move_grid %>%
          ungroup() %>%
          dplyr::select(from, to, prob_move) %>%
          # group_by(to) %>%
          # mutate(grouped_id = row_number()) %>%
          # ungroup()%>%
          spread(to, prob_move) %>%
          dplyr:: select(-from) %>%
          as.matrix()
        
      }
      
      if (num_patches > 1) {
        pop[now_year &
              pop$age > fish$min_age, ] <-
          move_fish(
            pop %>% filter(year == y, age > fish$min_age),
            fish = fish,
            num_patches = num_patches,
            move_matrix = adult_move_matrix
          )
      }
      
      # change management
      
      if ((y - burn_years) == manager$year_mpa) {
        pop$mpa[pop$patch %in% mpa_locations & pop$year >= y] <- T
        
        if (fleet$mpa_reaction == "leave") {
          mpa_effort <-
            sum(pop$effort[pop$patch %in% mpa_locations &
                             pop$year == (y - 1) & pop$age == 0])
          
          effort[y - 1] <-   effort[y - 1] - mpa_effort
          
          
        }
        
      }
      # Adjust fleet
      if (y > (burn_years)) {
        
        b0 <- pop %>%
          filter(year %in% seq(burn_years - 10, burn_years, 1)) %>%
          group_by(year) %>%
          summarise(b = sum(biomass),.groups="drop") %>%
          #ungroup() %>%
          summarise(b0 = mean(b))
        
        b0 <- b0$b0
        
        
        if (fleet$fleet_model == "open-access" &
            tune_costs == TRUE) {
          #located here to allow for b0
          
          # huh <- map_dbl(seq(0.01,0.9, by = 0.1), estimate_costs, fish = fish,
          #                fleet = fleet,
          #                b_ref_oa = fleet$b_ref_oa,
          #                lags = 10,
          #                num_patches = num_patches,
          #                sim_years = sim_years,
          #                burn_years = burn_years
          # )
          
          tuned_cr_ratio <-
            nlminb(
              c(0.3),
              estimate_costs,
              fish = fish,
              fleet = fleet,
              b_ref_oa = fleet$b_ref_oa,
              lower = c(1e-3,1e-3),
              upper = c(0.75,10),
              lags = fleet$profit_lags,
              num_patches = num_patches,
              sim_years = sim_years,
              burn_years = burn_years,
              sprinkler = sprinkler,
              mpa_habfactor = mpa_habfactor
            )
          
          # tuned_cr_ratio <-
          #   nlminb(
          #     c(0.3,0.05),
          #     estimate_costs,
          #     fish = fish,
          #     fleet = fleet,
          #     b_ref_oa = fleet$b_ref_oa,
          #     lower = c(1e-3,1e-3),
          #     upper = c(0.75,10),
          #     lags = 10,
          #     num_patches = num_patches,
          #     sim_years = sim_years,
          #     burn_years = burn_years,
          #     sprinkler = sprinkler,
          #     mpa_habfactor = mpa_habfactor
          #   )
          tune_costs <- FALSE
          
          fleet$max_cr_ratio <- tuned_cr_ratio$par[1]
          
          # fleet$max_perc_change_f <- tuned_cr_ratio$par[2]
          
          hyp_f <- fish$m #hypothetical f
          
          hyp_effort <- hyp_f / mean(q[(burn_years + 1):sim_years])
          
          hyp_f_at_age <- hyp_f * fleet$sel_at_age
          
          hyp_b0_catch <-
            sum((hyp_f_at_age / (hyp_f_at_age + fish$m))  * b0_at_age * (1 - exp(-(
              hyp_f_at_age + fish$m
            ))))
          
          b0_revenue <-
            mean(price_series[(burn_years + 1):sim_years]) * hyp_b0_catch
          
          hyp_profits_guess <- b0_revenue * (1 - fleet$max_cr_ratio)
          
          cost_guess <-
            (b0_revenue - hyp_profits_guess) / hyp_effort ^ fleet$beta
          
          fleet$theta <-
            (fleet$max_perc_change_f * hyp_effort) / (hyp_profits_guess / hyp_effort)
          
          cost_series <-
            generate_timeseries(
              cost_guess,
              cv = fleet$cost_cv,
              ac = fleet$cost_ac,
              percent_slope = fleet$cost_slope,
              time = sim_years
            )
          
          cost_series <- (cost_series / max(cost_series)) * cost_guess
          
          fleet$cost <- cost_guess
          
          cost_frame <- dplyr::data_frame(year = 1:sim_years, cost = cost_series)
          
          pop <- pop %>%
            dplyr::select(-cost) %>%
            dplyr::left_join(cost_frame, by = "year")
          
          if (fleet$cost_function == "distance from port") {
            cost_frame <-
              expand.grid(year = 1:sim_years, patch = 1:num_patches) %>%
              dplyr::as_data_frame() %>%
              dplyr::left_join(cost_frame, by = "year") %>%
              dplyr::mutate(cost = cost * (1 + fleet$cost_slope * (patch - 1)))
            
            pop <- pop %>%
              dplyr::select(-cost) %>%
              dplyr::left_join(cost_frame, by = c("patch", "year"))
          }
          
        }
        
        if (y == (burn_years + 2) &
            fleet$fleet_model == "open-access") {
          profits <- pop %>%
            filter(year >= (y - (1 + fleet$profit_lags)), year < y) %>%
            group_by(year) %>%
            summarise(profits = sum(profits),.groups="drop")
          
          total_initial_profits <- mean(profits$profits)
          
        }
        
        previous_max <-
          ifelse(y > (burn_years + 1), max(effort[max(1, (y - 1 - max_window)):(y - 1)]), fleet$initial_effort)
        
        effort[y] <- determine_effort(
            last_effort = ifelse(y > (burn_years + 1), effort[y - 1], fleet$initial_effort),
            fleet = fleet,
            fish = fish,
            y = y,
            burn_years = burn_years,
            pop = pop,
            mpa = mpa,
            num_patches = num_patches,
            effort_devs = effort_devs,
            profit_lags = fleet$profit_lags,
            previous_max = previous_max
        )
        
      }
      
      pop[now_year, "effort"] <-
        distribute_fleet(
          pop = pop %>% filter(year == y),
          prior_profits = pop$profits[pop$year == (y - 1)],
          prior_effort = pop$effort[pop$year == (y - 1)],
          year = y,
          burn_years = burn_years,
          effort = effort[y],
          fleet = fleet,
          num_patches = num_patches,
          mpa = mpa
        )
      
      if (fleet$tech_rate > 0 & y > burn_years) {
        q[y] <-
          q[y - 1] + pmax(0,
                          rnorm(1, fleet$tech_rate * q[y - 1], fleet$tech_rate * fleet$q))
      }
      
      pop[now_year, "f"] <-
        pop[now_year, "effort"] * q[y]
      
      # grow and die -----
      
      if(keep_burn==TRUE) {
        kap <- sum(pop[pop$year==burn_years,"numbers"]) ## based on abundance, not biomass
      }
      
      pop[pop$year == (y + 1), "numbers"] <-
        pop[now_year, ] %>%
        group_by(patch) %>%
        mutate(totPop=sum(biomass))%>% # calc total pop for log growth calcs
        dplyr::mutate(numbers = 
                        ifelse(patch %in% mpa_locations,
                               grow_and_die_adj(
                                 numbers = numbers,
                                 f = f,
                                 mpa = TRUE,
                                 fish = fish,
                                 fleet = fleet,
                                 y = y,
                                 mImp=mImp,
                                 farmYrs=farmYrs)$survivors,
                               grow_and_die_adj(
                                 numbers = numbers,
                                 f = f,
                                 mpa = FALSE,
                                 fish = fish,
                                 fleet = fleet,
                                 mImp=1,
                                 farmYrs=farmYrs)$survivors)) %>%
        ungroup() %>%
        {
          .$numbers
        }
      
      pop[now_year, "numbers_caught"] <-
        pop[now_year, ] %>%
        group_by(patch) %>%
        mutate(totPop=sum(biomass))%>% # calc total pop for log growth calcs
        dplyr::mutate(
          numbers_caught =
            ifelse(patch %in% mpa_locations,
                   c(grow_and_die_adj(
                     numbers = numbers,
                     f = f,
                     mpa = TRUE,
                     fish = fish,
                     fleet = fleet,
                     y = y,
                     mImp=mImp,
                     farmYrs=farmYrs)$caught),
                   c(grow_and_die_adj(
                     numbers = numbers,
                     f = f,
                     mpa = FALSE,
                     fish = fish,
                     fleet = fleet,
                     y = y,
                     mImp=1,
                     farmYrs=farmYrs)$caught))) %>%
        ungroup() %>%
        {
          .$numbers_caught
        }
      
      pBmDf <- pop%>%
        filter(year==y)%>%
        left_join(.,pop%>%
                    filter(year == y-1)%>%
                    group_by(patch)%>%
                    summarise(pBm=sum(biomass),.groups="drop"),
                  by="patch")%>%
        mutate(biomass2=as.numeric(purrr::pmap(list(age,numbers,patch,pBm,mpa),~waak(aGe=..1,
                                                                                nums=..2,
                                                                                ptch=..3, 
                                                                                mpa = ..5,
                                                                                y=y,
                                                                                farmYrs=farmYrs,
                                                                                pI=prodIncr,
                                                                                pBm =..4,
                                                                                pK=patchK,
                                                                                ddwaa=ddwaa))),
               biomass_caught2 = as.numeric(purrr::pmap(list(age,numbers_caught,patch,pBm,mpa),~waak(aGe=..1,
                                                                                         nums=..2,
                                                                                         ptch=..3, 
                                                                                         mpa = ..5,
                                                                                         y=y,
                                                                                         farmYrs=farmYrs,
                                                                                         pI=prodIncr,
                                                                                         pBm =..4,
                                                                                         pK=patchK,
                                                                                         ddwaa=ddwaa))))%>%
        select(year,age,patch,biomass2,biomass_caught2)

     
        pop<-pop%>%
          left_join(.,pBmDf,by=c("year","age","patch"))%>%
          mutate(biomass=case_when(is.na(biomass2) ~ biomass,
                                   TRUE ~ biomass2),
                 biomass_caught=case_when(is.na(biomass_caught2) ~ biomass_caught,
                                   TRUE ~ biomass_caught2))

      pop <- pop %>%
        dplyr::mutate(patch_age_costs = ((cost) * (effort) ^ fleet$beta) / fish$max_age) %>% # divide costs up among each age class
        dplyr::mutate(
          ssb = numbers * ssb_at_age,
          #numbers * weight_at_age,
          # biomass_caught = numbers_caught * weight_at_age,
          profits = biomass_caught * price - patch_age_costs
        )
      
      # spawn ----

            # if (is.na(spawning_season) | ((((year) - floor(year))/spawning_season) == 1))
      pop$numbers[pop$year == (y + 1) &
                    pop$age == fish$min_age] <-
        calculate_recruits(
          #calculate_recruitsFarm(
          #farmYrs=farmYrs,
          pop = pop[pop$year == y, ],
          fish = fish,
          num_patches = num_patches,
          phase = model_phase,
          move_matrix = larval_move_matrix,
          rec_devs = rec_devs[y + 1],
          patch_habitat = habitat
        )
      
      
      if (y == burn_years) {
        fish$ssb0 <- pop %>%
          filter(year == burn_years) %>%
          group_by(patch) %>%
          summarise(ssb = sum(ssb),.groups="drop") %>%
          ungroup() %>%
          {
            (.$ssb)
          }
        
        model_phase <- "recruit"
        
        effort[y + 1] <- fleet$initial_effort
      }
    }
    rec_mat <- dplyr::data_frame(year = 1:sim_years, rec_dev = rec_devs)
    
    enviro_mat <- dplyr::data_frame(year = 1:sim_years, enviro = enviro)
    
    og <- burn_years
    if (keep_burn == TRUE) {
      burn_years <- -99
    }
    
    pop <- pop %>%
      dplyr::left_join(rec_mat, by = "year") %>%
      dplyr::left_join(enviro_mat, by = "year") %>%
      dplyr::filter(year > burn_years, year < max(year)) %>%
      dplyr:: mutate(
        burn = year <= og,
        eventual_mpa = patch %in% mpa_locations,
        msy = msy,
        e_msy = e_msy,
        b0 = b0
      )
    
    return(pop)
  }


#####################################################################
#####################################################################

######### grow and die adjusted: adds mortality impacts (mImp) and carrying capacity growth, which req pK input (as 'kap' in simFisheryAqK)

grow_and_die_adj <- function(numbers, f, mpa, fish, fleet, y, mImp,farmYrs=farmYrs) {
  
  if(mpa == TRUE & farmYrs[y] == TRUE) {

    mImp <- mImp
  } else { 
    mImp <-1

    }
  
  survivors <- vector(mode = 'numeric', length = length(numbers))
  
  survival<-(exp(-fish$time_step*(fish$m*(mImp) + (f * fleet$sel_at_age))))

    death <-  1 - survival
  
  max_index <- length(survivors)


  survivors[2:max_index] <-
    numbers[1:(max_index - 1)] * survival[1:(max_index - 1)]

  survivors[max_index] <-
    survivors[max_index] +  numbers[max_index] * survival[max_index]

  
  caught <-((fish$time_step * f * fleet$sel_at_age) / (fish$time_step * (fish$m*(mImp) + (f * fleet$sel_at_age)))) *  (numbers * death)
  
  
  return(list(survivors = survivors, caught = caught))

}

####################################################################################
####################################################################################

waak <- function(aGe, nums, ptch, mpa, farmYrs=farmYrs, y=y, pBm=pBm, pK=patchK, pI=prodIncr,ddwaa=ddwaa){

  if(mpa == TRUE & farmYrs[y] == TRUE) {
    pK <- pI*pK
  } else { 
    pK <- pK
  }

  pGr <- ddwaa$statGrwth[ddwaa$age == aGe] * (1-(pBm/pK)) 

  pWt <- if(aGe==1) { 0 } else {ddwaa$waa[ddwaa$age == aGe-1] + pGr}
  
  bioM <- nums * pWt
  
  return(bioM) #biomass by age & patch
}
