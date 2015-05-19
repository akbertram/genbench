WITH good_ids AS (
            select 
              pubmed_ids 
            from 
              generifs_basic
            where 
              tax_id=9606
              AND NOT annotation LIKE '%(HuGE Navigator)' -- consider deleting for full RIFome
            group by 
              pubmed_ids, annotation
            having count(distinct gene_id) > 1

)

select distinct 
  rif.pubmed_ids,
  rif.annotation,
  rif.gene_id
from
  generifs_basic rif
where
 pubmed_ids in good_ids
order by pubmed_ids;


