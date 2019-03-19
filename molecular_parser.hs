import Data.List
import Data.Array

bondDict n | n > 1 = ["=#$" !! (n - 2)]
           | otherwise = ""

aceticAcid = ([[0, 1, 2, 1],
               [1, 0, 0, 0],
               [2, 0, 0, 0],
               [1, 0, 0, 0]], ["C", "C", "O", "O"])

aceticAcidSMILES = "CC(O)=O"

cyclicMol = ([[0, 1, 0, 0, 1, 0],
              [1, 0, 1, 0, 0, 0],
              [0, 1, 0, 1, 0, 0],
              [0, 0, 1, 0, 1, 0],
              [1, 0, 0, 1, 0, 1],
              [0, 0, 0, 0, 1, 0]] , ["C", "C", "C", "N", "C", "O"])

cmEL = adjToEdges (fst cyclicMol)
aaEL = adjToEdges (fst aceticAcid)

simpleEdgeList [] = []
simpleEdgeList (e:es) = if elem e es then simpleEdgeList es else e: simpleEdgeList es

edgeListToAdjList [] = []
edgeListToAdjList el = zip [0..m] [[snd e | e <- el, fst e == i] | i <- [0..m]]
    where m = maximum $ concat [[fst x, snd x] | x <- el]

adjToEdges []     = []
adjToEdges rs = convert 0 0 rs

convert n k [] = []
convert n k (r:rs) = concat [multiEdge c (n, ind) | (c, ind) <- zip r [0..], c > k, ind > n]
                     ++ convert (n+1) k rs

multiEdge 0 e = []
multiEdge n e =  e:(multiEdge (n-1) e)

-- build up a minmum spanning tree
--prims :: [(Int, Int)] -> [(Int, Int)]
prims [] = []
prims el = sort $ prims' m [0] [] el
      where m = maximum $ concat [[fst x, snd x] | x <- el]

prims' m ns es el | genericLength ns == m + 1 = es
                  | otherwise = prims' m ((fst e):ns) ((snd e):es) el
                      where e = findEdge ns el

--findEdgePrim      :: [Int] -> [(Int, Int)] -> (Int, (Int, Int))
findEdge ns [] = error "No edges found"
findEdge ns (e:es) | (elem (fst e) ns) && (notElem (snd e) ns) = ((snd e), e)
                   | (notElem (fst e) ns) && (elem (snd e) ns) = ((fst e), e)
                   | otherwise = findEdge ns es

cyclicEdges el = [e | e <- el, notElem e p]
      where p = prims el

bondCount el = [countBond e el | e <- p]
      where p = prims el

countBond e el = (c, e)
      where c = genericLength [e' | e' <- el, e' == e]

contents i []                  = []
contents i (a:as) | i == fst a = snd a
                  | otherwise  = contents i as

--cycleFinder :: [(Int, Int)] -> [(Int, Int)] -> [[(Int)]]
cycleFinder [] el = []
cycleFinder (c:cs) el = (findSmallestCycle' [fst c] (snd c) [c] el') : cycleFinder cs el
      where el' = removeEdge c el

--findSmallestCycle' :: Ord a => [a] -> a -> [(a, a)] -> [(a, a)] -> [(a, a)]
findSmallestCycle' cv end ce el = head . sort $ zip (map genericLength z) z
      where z = findSmallestCycle cv end ce el

--findSmallestCycle :: Int -> Int -> [(Int, Int)] -> [(Int)]
-- end current_edges remaining_edge_list
findSmallestCycle cv end ce [] = [ce]
findSmallestCycle cv end ce el = concat [continueCycle cv' end ce' el' |
                                  n <- nodes,
                                  let cv' = nub $ (fst n):(snd n):cv,
                                  let ce' = n:ce,
                                  let el' = removeEdge n el]
      where nodes = nub $ [e | e <- el, elem (fst e) cv || elem (snd e) cv]

continueCycle cv end ce el | elem end cv = [ce]
                           | otherwise   = findSmallestCycle cv end ce el

removeEdge e [] = []
removeEdge e' (e:es) = if e == e' then es else e : removeEdge e' es


-- Graph to String Driver Functions
makeSMILES mol = smilesify mol 0

smilesify mol n = key !! n ++ ringChk n cyc bc ++ smileMore n (contents n al) mol
      where key = snd mol
            el = adjToEdges (fst mol)
            el' = simpleEdgeList el
            p = prims el
            al = edgeListToAdjList p
            bc  = bondCount el
            cyc = cyclicEdges el

smileMore p [] mol     = ""
smileMore p (v:vs) mol | null vs   = bondChk (p, v) bc ++ smilesify mol v
                       | otherwise =  "(" ++ bondChk (p, v) bc ++ smilesify mol v ++ ")" ++ smileMore p vs mol
                        where el = adjToEdges (fst mol)
                              bc  = bondCount el
                              cyc = cyclicEdges el

bondChk e []                  = ""
bondChk e (b:bs) | e == snd b = bondDict (fst b)
                 | otherwise  = bondChk e bs

ringChk n [] _               = ""
ringChk n cyc bc | elem n v  = concat
                              [show i ++ bondChk e bc |
                              (i, e) <- zip [0..] cyc, n == fst e || n == snd e]
                 | otherwise = ""
                 where v = nub (map fst cyc ++ map snd cyc)

--huckels :: Num a => [String] -> [(a, a)] -> Bool
huckels key [] bc = False
huckels key cs bc = (length [i | i <- [1..c], 4 * i + 2 == c]) > 1
    where c = huckelAtomCount key cs + huckelBondCount bc

huckelAtomCount [] key = 0
huckelAtomCount (c:cs) key | elem (key !! c) ["O", "N", "S"]  = 1 + huckelAtomCount cs key
                           | otherwise                        = 0 + huckelAtomCount cs key

huckelBondCount [] = 0
huckelBondCount (b:bs) | b == "="  = 1 + huckelBondCount bs
                       | b == "#"  = 2 + huckelBondCount bs
                       | b == "$"  = 3 + huckelBondCount bs
                       | otherwise = 0
