

-- Query to match gene ids with assay ids
select * from intended_target IT, gene G where G.gene_id=IT.target_id;


-- Match gsid and sample id in sample table
-- To be used in any query
select sa_gsid, sa_sample_id from sample 
where sa_gsid = '';

-- Match aeid to assay name 
select * from assay_component_endpoint
where  aeid = '';


-- Extract hit calls for a specific spid 
select spid ,l5.aeid, hitc, modl_ga from level5 l5, level4 l4
where l5.l4id = l4.l4id 
and l4.spid = '' ;
