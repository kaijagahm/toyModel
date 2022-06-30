makeGraphsList<- function (dataList, weighted = FALSE, id1Col = "ID1", id2Col = "ID2") 
{
  simplified <- lapply(dataList, function(x) {
    checkmate::assertChoice("interval", names(x))
    x <- x %>% dplyr::select(.data[[id1Col]], .data[[id2Col]], interval)
  })
  if (weighted == FALSE) {
    simplified <- lapply(simplified, function(x) {
      x <- x %>% dplyr::distinct()
    })
    gs <- lapply(simplified, function(x) {
      igraph::graph_from_data_frame(d = x, directed = FALSE)
    })
    names(gs) <- unlist(lapply(simplified, function(x){unique(x$interval)}))
  }
  else {
    simplified <- lapply(simplified, function(x) {
      x <- x %>% dplyr::mutate(weight = 1) %>% dplyr::group_by(.data[[id1Col]], 
                                                               .data[[id2Col]]) %>% dplyr::summarize(weight = sum(weight)) %>% 
        dplyr::ungroup()
    })
    gs <- lapply(simplified, function(x) {
      igraph::graph_from_data_frame(d = x, directed = FALSE)
    })
    names(gs) <- unlist(lapply(simplified, function(x){unique(x$interval)}))
    
  }
  return(list(graphs = gs, simplifiedData = simplified))
}




makeGraphs <- function (edges, fullData, interval, dateTimeStart = NULL, dateTimeEnd = NULL, 
                        id1Col = "ID1", id2Col = "ID2", weighted = FALSE) 
{
  checkmate::assertLogical(weighted, len = 1)
  checkmate::assertDataFrame(fullData)
  checkmate::assertChoice("timegroup", names(fullData))
  checkmate::assertChoice("timestamp", names(fullData))
  timegroupInfo <- fullData %>% dplyr::select(timegroup, timestamp) %>% 
    dplyr::group_by(timegroup) %>% dplyr::summarize(minDatetime = min(timestamp), 
                                                    maxDatetime = max(timestamp))
  int <- lubridate::as.duration(interval)
  if (is.na(int)) {
    stop("Argument `interval` could not be expressed as a duration: lubridate::as.duration() returned NA. Please make sure you are specifying a valid interval, such as '1 day', '3 hours', '2 weeks', etc.")
  }
  checkmate::assertClass(int, "Duration")
  if (is.null(dateTimeStart)) {
    dateTimeStart <- min(fullData$timestamp)
    warning(paste0("No start datetime provided. Using earliest `timestamp` from `fullData`, which is ", 
                   dateTimeStart, "."))
  }
  if (is.null(dateTimeEnd)) {
    dateTimeEnd <- max(fullData$timestamp)
    warning(paste0("No end datetime provided. Using latest `timestamp` from `fullData`, which is ", 
                   dateTimeEnd, "."))
  }
  start <- lubridate::parse_date_time(dateTimeStart, orders = c("%Y%m%d %H%M%S", 
                                                                "%Y%m%d %H%M", "%Y%m%d"))
  if (is.na(start)) {
    stop("`dateTimeStart` could not be parsed. Please make sure you have used one of the following formats: YYYY-MM-DD hh:mm:ss, YYYY-MM-DD hh:mm, or YYYY-MM-DD.")
  }
  end <- lubridate::parse_date_time(dateTimeEnd, orders = c("%Y%m%d %H%M%S", 
                                                            "%Y%m%d %H%M", "%Y%m%d"))
  if (is.na(end)) {
    stop("`dateTimeEnd` could not be parsed. Please make sure you have used one of the following formats: YYYY-MM-DD hh:mm:ss, YYYY-MM-DD hh:mm, or YYYY-MM-DD.")
  }
  timegroupInfo <- timegroupInfo %>% tibble::add_row(minDatetime = start, 
                                                     .before = 1) %>% tibble::add_row(minDatetime = end)
  breaks <- seq(from = start, to = end, by = int)
  groupedTimegroups <- timegroupInfo %>% dplyr::mutate(interval = cut(minDatetime, 
                                                                      breaks)) %>% dplyr::select(timegroup, interval)
  checkmate::assertDataFrame(edges)
  checkmate::assertChoice(id1Col, names(edges))
  checkmate::assertChoice(id2Col, names(edges))
  checkmate::assertChoice("timegroup", names(edges))
  dataList <- edges %>% ungroup() %>% dplyr::select(.data[[id1Col]], 
                                                    .data[[id2Col]], timegroup) %>% dplyr::left_join(groupedTimegroups, 
                                                                                                     by = "timegroup") %>% dplyr::group_by(interval) %>% 
    dplyr::group_split(.keep = TRUE)

  networks <- makeGraphsList(dataList = dataList, 
                                           weighted = weighted, id1Col = id1Col, id2Col = id2Col)
  return(networks)
}
