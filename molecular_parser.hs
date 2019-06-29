import Data.List
import Data.Array
import Data.Char

bondDict n | n > 1 = ["=#$" !! (n - 2)]
           | otherwise = ""

aceticAcid = ([[0, 1, 2, 1],
               [1, 0, 0, 0],
               [2, 0, 0, 0],
               [1, 0, 0, 0]], ["C", "C", "O", "O"])

aceticAcidSMILES = "CC(O)=O"

cyclopentanol = ([[0, 1, 0, 0, 1, 0],
                  [1, 0, 1, 0, 0, 0],
                  [0, 1, 0, 1, 0, 0],
                  [0, 0, 1, 0, 1, 0],
                  [1, 0, 0, 1, 0, 1],
                  [0, 0, 0, 0, 1, 0]] ,
                  ["C", "C", "C", "C", "C", "O"])

benzene = ([[0,2,0,0,0,1],
            [2,0,1,0,0,0],
            [0,1,0,2,0,0],
            [0,0,2,0,1,0],
            [0,0,0,1,0,2],
            [2,0,0,0,1,0]], ["C", "C", "C", "C", "C", "C"])

cycloEdges = adjToEdges (fst cyclopentanol)
aceticEdges = adjToEdges (fst aceticAcid)

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

spanningTree [] = []
spanningTree el = sort $ spanningTree' m [0] [] el
      where m = maximum $ mergeTupleList el

spanningTree' m ns es el | genericLength ns == m + 1 = es
                         | otherwise = spanningTree' m ((fst e):ns) ((snd e):es) el
      where e = findEdge ns el

findEdge ns [] = error "No edges found"
findEdge ns (e:es) | (elem (fst e) ns) && (notElem (snd e) ns) = ((snd e), e)
                   | (notElem (fst e) ns) && (elem (snd e) ns) = ((fst e), e)
                   | otherwise = findEdge ns es

cyclicEdges el = [e | e <- el, notElem e st]
      where st = spanningTree el

bondCount el = [countBond e el | e <- st]
      where st = spanningTree el

countBond e el = (c, e)
      where c = genericLength (filter (== e) el)

cycleFinder [] el = []
cycleFinder (c:cs) el = (findSmallestCycle' [fst c] (snd c) [c] el') : cycleFinder cs el
      where el' = removeEdge c el

findSmallestCycle' cv end ce el = head . sort $ zip (map genericLength z) z
      where z = findSmallestCycle cv end ce el

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

removeEdge e el = filter (/= e) el

-- Graph to String Driver Functions
makeSMILES mol = smilesify mol 0

smilesify mol n = key !! n ++ ringChk n cyc bc ++ smileMore n (contents n al) mol
      where key = snd mol
            el = adjToEdges (fst mol)
            el' = simpleEdgeList el
            st = spanningTree el
            al = edgeListToAdjList st
            bc  = bondCount el
            cyc = cyclicEdges el'

smileMore st [] mol     = ""
smileMore st (v:vs) mol | null vs   = bondChk (st, v) bc ++ smilesify mol v
                        | otherwise = "(" ++ bondChk (st, v) bc ++ smilesify mol v ++ ")"
                                     ++ smileMore st vs mol
                        where el = adjToEdges (fst mol)
                              bc  = bondCount el
                              cyc = cyclicEdges el

bondChk e []                  = ""
bondChk e (b:bs) | e == snd b = bondDict (fst b)
                 | otherwise  = bondChk e bs

ringChk n [] _               = ""
ringChk n cyc bc | elem n v  = concat
                              [show i ++ bondChk e bc |
                              (i, e) <- zip [1..] cyc, n == fst e || n == snd e]
                 | otherwise = ""
                 where v = mergeTupleList cyc

contents i []                  = []
contents i (a:as) | i == fst a = snd a
                 | otherwise  = contents i as

mergeTupleList xs      = nub $ mergeTupleList' xs
mergeTupleList' []     = []
mergeTupleList' (x:xs) = (fst x):(snd x):(mergeTupleList xs)

{-

-- Nonfunctioning aromatic conversion code, typing issues

convertAromatics (am, key) = [if elem i arm
                              then (map toLower k)
                              else k | (i, k) <- zip [0..] key]
                              where el  = adjToEdges am
                                    el' = simpleEdgeList el
                                    bc  = bondCount el
                                    cycE = cyclicEdges el'
                                    cycles = cycleFinder cycE el'
                                    arm = concat [aromatics (snd c) bc key | c <- cycles]

-- cyclic edge tuple + bond counts + key string
-- checks if this ring is aromatic returns a list of aromatic vertices to conver if true
aromatics cyc bc key = if conjugated cyc bc key
                       then
                          if huckels key cyc bc
                            then mergeTupleList cyc
                            else []
                          else []

conjugated cyc bc key = conjugated' vl (snd cyc) bc key
    where vl = mergeTupleList cyc

conjugated' []     cyc bc key = True -- all atoms in cycle are conjugated
conjugated' (v:vl) cyc bc key | elem (key !! v) ["O", "N", "P", "S"]  = conjugated' vl cyc bc key -- atom is a heteroatom
                              | elem v doubleBondVerts                = conjugated' vl cyc bc key -- atom is part of a double bond
                              | otherwise                             = False   -- structure not conjugated
                              where doubleBondVerts = mergeTupleList [(snd a) | a <- bc, (fst a) > 1]


huckels key [] bc = False
huckels key cs bc = (genericLength [i | i <- [1..c], 4 * i + 2 == c]) >= 1
    where c = huckelAtomCount key cs + huckelBondCount bc

huckelAtomCount [] key = 0
huckelAtomCount (c:cs) key | elem (key !! c) ["O", "N", "S", "P"]  = 1 + huckelAtomCount cs key
                           | otherwise                             = huckelAtomCount cs key

huckelBondCount [] = 0
huckelBondCount (b:bs) | fst b == 2  = 2 + huckelBondCount bs -- double
                       | fst b == 3  = 4 + huckelBondCount bs -- triple
                       | fst b == 4  = 4 + huckelBondCount bs -- quadruple
                       | otherwise = 0
-}
