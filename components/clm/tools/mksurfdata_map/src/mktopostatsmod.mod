
  %  R   k820309    a          16.0        :��X                                                                                                           
       mktopostatsMod.F90 MKTOPOSTATSMOD              MKTOPOSTATS                      @                              
       SHR_SYS_FLUSH                      @                              
       DOMAIN_CHECKSAME                                                     
       R8 SHR_KIND_R8                   @                               '�                   #SET    #NS    #NI    #NJ    #EDGEN 	   #EDGEE 
   #EDGES    #EDGEW    #MASK    #FRAC    #LATC    #LONC    #LATS    #LATN    #LONW    #LONE    #AREA    #IS_2D    #FRACSET    #MASKSET                 � $                                                                      � $                                                             � $                                                             � $                                                             � $                             	                
                � $                             
     (          
                � $                                  0          
                � $                                  8          
               �$                                          @              	               &                                                       �$                                         �              
   
            &                                                       �$                                         �                 
            &                                                       �$                                                         
            &                                                       �$                                         `                
            &                                                       �$                                         �                
            &                                                       �$                                         �                
            &                                                       �$                                         8                
            &                                                       �$                                         �                
            &                                                        � $                                   �                         � $                                   �                         � $                                   �                              @                               '                    #SET    #NAME    #NA    #NB    #NI    #NJ    #NS     #YC_SRC !   #YC_DST "   #XC_SRC #   #XC_DST $   #MASK_SRC %   #MASK_DST &   #AREA_SRC '   #AREA_DST (   #FRAC_SRC )   #FRAC_DST *   #SRC_INDX +   #DST_INDX ,   #WOVR -                � $                                                                      � $                                                                      � $                                   @                          � $                                   D                          � $                                   H                          � $                                   L                          � $                                    P                         �$                             !            X                 
            &                                                       �$                             "            �              	   
            &                                                       �$                             #            �              
   
            &                                                       �$                             $            0                
            &                                                       �$                              %            x                            &                                                       �$                              &            �                            &                                                       �$                             '                            
            &                                                       �$                             (            P                
            &                                                       �$                             )            �                
            &                                                       �$                             *            �                
            &                                                       �$                              +            (                            &                                                       �$                              ,            p                            &                                                       �$                             -            �                
            &                                                                                        .                                                         #         @                                  /                    #UNIT 0                                             0            #         @                                  1                    #SRCDOMAIN 2   #DSTDOMAIN 3   #TGRIDMAP 4             
                                  2     �             #DOMAIN_TYPE              
                                  3     �             #DOMAIN_TYPE              
                                  4                   #GRIDMAP_TYPE    #         @                                   5                   #MKTOPOSTATS%DOMAIN_TYPE 6   #LDOMAIN K   #MAPFNAME L   #DATFNAME M   #NDIAG N   #TOPO_STDDEV_O O   #SLOPE_O P                     @                           6     '�                   #SET 7   #NS 8   #NI 9   #NJ :   #EDGEN ;   #EDGEE <   #EDGES =   #EDGEW >   #MASK ?   #FRAC @   #LATC A   #LONC B   #LATS C   #LATN D   #LONW E   #LONE F   #AREA G   #IS_2D H   #FRACSET I   #MASKSET J                � $                              7                                        � $                              8                               � $                              9                               � $                              :                               � $                             ;                
                � $                             <     (          
                � $                             =     0          
                � $                             >     8          
               �$                              ?            @              	               &                                                       �$                             @            �              
   
            &                                                       �$                             A            �                 
            &                                                       �$                             B                            
            &                                                       �$                             C            `                
            &                                                       �$                             D            �                
            &                                                       �$                             E            �                
            &                                                       �$                             F            8                
            &                                                       �$                             G            �                
            &                                                        � $                              H     �                         � $                              I     �                         � $                              J     �                      
  @                               K     �             #MKTOPOSTATS%DOMAIN_TYPE 6             
  @                              L                    1           
@ @                              M                    1           
  @                               N                     D @                              O                   
               &                                                     D @                              P                   
               &                                              �   *      fn#fn $   �      b   uapp(MKTOPOSTATSMOD    �   N   J  SHR_SYS_MOD    4  Q   J  MKDOMAINMOD    �  O   J  SHR_KIND_MOD (   �        DOMAIN_TYPE+MKDOMAINMOD ,   �  P   a   DOMAIN_TYPE%SET+MKDOMAINMOD +   @  H   a   DOMAIN_TYPE%NS+MKDOMAINMOD +   �  H   a   DOMAIN_TYPE%NI+MKDOMAINMOD +   �  H   a   DOMAIN_TYPE%NJ+MKDOMAINMOD .     H   a   DOMAIN_TYPE%EDGEN+MKDOMAINMOD .   `  H   a   DOMAIN_TYPE%EDGEE+MKDOMAINMOD .   �  H   a   DOMAIN_TYPE%EDGES+MKDOMAINMOD .   �  H   a   DOMAIN_TYPE%EDGEW+MKDOMAINMOD -   8  �   a   DOMAIN_TYPE%MASK+MKDOMAINMOD -   �  �   a   DOMAIN_TYPE%FRAC+MKDOMAINMOD -   `  �   a   DOMAIN_TYPE%LATC+MKDOMAINMOD -   �  �   a   DOMAIN_TYPE%LONC+MKDOMAINMOD -   �  �   a   DOMAIN_TYPE%LATS+MKDOMAINMOD -     �   a   DOMAIN_TYPE%LATN+MKDOMAINMOD -   �  �   a   DOMAIN_TYPE%LONW+MKDOMAINMOD -   D	  �   a   DOMAIN_TYPE%LONE+MKDOMAINMOD -   �	  �   a   DOMAIN_TYPE%AREA+MKDOMAINMOD .   l
  H   a   DOMAIN_TYPE%IS_2D+MKDOMAINMOD 0   �
  H   a   DOMAIN_TYPE%FRACSET+MKDOMAINMOD 0   �
  H   a   DOMAIN_TYPE%MASKSET+MKDOMAINMOD *   D  5      GRIDMAP_TYPE+MKGRIDMAPMOD .   y  P   a   GRIDMAP_TYPE%SET+MKGRIDMAPMOD /   �  P   a   GRIDMAP_TYPE%NAME+MKGRIDMAPMOD -     H   a   GRIDMAP_TYPE%NA+MKGRIDMAPMOD -   a  H   a   GRIDMAP_TYPE%NB+MKGRIDMAPMOD -   �  H   a   GRIDMAP_TYPE%NI+MKGRIDMAPMOD -   �  H   a   GRIDMAP_TYPE%NJ+MKGRIDMAPMOD -   9  H   a   GRIDMAP_TYPE%NS+MKGRIDMAPMOD 1   �  �   a   GRIDMAP_TYPE%YC_SRC+MKGRIDMAPMOD 1     �   a   GRIDMAP_TYPE%YC_DST+MKGRIDMAPMOD 1   �  �   a   GRIDMAP_TYPE%XC_SRC+MKGRIDMAPMOD 1   =  �   a   GRIDMAP_TYPE%XC_DST+MKGRIDMAPMOD 3   �  �   a   GRIDMAP_TYPE%MASK_SRC+MKGRIDMAPMOD 3   e  �   a   GRIDMAP_TYPE%MASK_DST+MKGRIDMAPMOD 3   �  �   a   GRIDMAP_TYPE%AREA_SRC+MKGRIDMAPMOD 3   �  �   a   GRIDMAP_TYPE%AREA_DST+MKGRIDMAPMOD 3   !  �   a   GRIDMAP_TYPE%FRAC_SRC+MKGRIDMAPMOD 3   �  �   a   GRIDMAP_TYPE%FRAC_DST+MKGRIDMAPMOD 3   I  �   a   GRIDMAP_TYPE%SRC_INDX+MKGRIDMAPMOD 3   �  �   a   GRIDMAP_TYPE%DST_INDX+MKGRIDMAPMOD /   q  �   a   GRIDMAP_TYPE%WOVR+MKGRIDMAPMOD ,     p       R8+SHR_KIND_MOD=SHR_KIND_R8 *   u  R       SHR_SYS_FLUSH+SHR_SYS_MOD /   �  @   a   SHR_SYS_FLUSH%UNIT+SHR_SYS_MOD -     t       DOMAIN_CHECKSAME+MKDOMAINMOD 7   {  Y   a   DOMAIN_CHECKSAME%SRCDOMAIN+MKDOMAINMOD 7   �  Y   a   DOMAIN_CHECKSAME%DSTDOMAIN+MKDOMAINMOD 6   -  Z   a   DOMAIN_CHECKSAME%TGRIDMAP+MKDOMAINMOD    �  �       MKTOPOSTATS 4   @       MKTOPOSTATS%DOMAIN_TYPE+MKDOMAINMOD 8   \  P   a   MKTOPOSTATS%DOMAIN_TYPE%SET+MKDOMAINMOD 7   �  H   a   MKTOPOSTATS%DOMAIN_TYPE%NS+MKDOMAINMOD 7   �  H   a   MKTOPOSTATS%DOMAIN_TYPE%NI+MKDOMAINMOD 7   <  H   a   MKTOPOSTATS%DOMAIN_TYPE%NJ+MKDOMAINMOD :   �  H   a   MKTOPOSTATS%DOMAIN_TYPE%EDGEN+MKDOMAINMOD :   �  H   a   MKTOPOSTATS%DOMAIN_TYPE%EDGEE+MKDOMAINMOD :     H   a   MKTOPOSTATS%DOMAIN_TYPE%EDGES+MKDOMAINMOD :   \  H   a   MKTOPOSTATS%DOMAIN_TYPE%EDGEW+MKDOMAINMOD 9   �  �   a   MKTOPOSTATS%DOMAIN_TYPE%MASK+MKDOMAINMOD 9   8  �   a   MKTOPOSTATS%DOMAIN_TYPE%FRAC+MKDOMAINMOD 9   �  �   a   MKTOPOSTATS%DOMAIN_TYPE%LATC+MKDOMAINMOD 9   `  �   a   MKTOPOSTATS%DOMAIN_TYPE%LONC+MKDOMAINMOD 9   �  �   a   MKTOPOSTATS%DOMAIN_TYPE%LATS+MKDOMAINMOD 9   �  �   a   MKTOPOSTATS%DOMAIN_TYPE%LATN+MKDOMAINMOD 9      �   a   MKTOPOSTATS%DOMAIN_TYPE%LONW+MKDOMAINMOD 9   �   �   a   MKTOPOSTATS%DOMAIN_TYPE%LONE+MKDOMAINMOD 9   D!  �   a   MKTOPOSTATS%DOMAIN_TYPE%AREA+MKDOMAINMOD :   �!  H   a   MKTOPOSTATS%DOMAIN_TYPE%IS_2D+MKDOMAINMOD <    "  H   a   MKTOPOSTATS%DOMAIN_TYPE%FRACSET+MKDOMAINMOD <   h"  H   a   MKTOPOSTATS%DOMAIN_TYPE%MASKSET+MKDOMAINMOD $   �"  e   a   MKTOPOSTATS%LDOMAIN %   #  L   a   MKTOPOSTATS%MAPFNAME %   a#  L   a   MKTOPOSTATS%DATFNAME "   �#  @   a   MKTOPOSTATS%NDIAG *   �#  �   a   MKTOPOSTATS%TOPO_STDDEV_O $   y$  �   a   MKTOPOSTATS%SLOPE_O 