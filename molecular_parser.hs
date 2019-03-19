import Data.Graph
import Data.List
import Data.Array


bondDict n | n > 1 = ["=#$" !! (n - 2)]
           | otherwise = ""

--aceticAcid :: ([[Int]], [Atom])
aceticAcid = ([[0, 1, 2, 1],
               [1, 0, 0, 0],
               [2, 0, 0, 0],
               [1, 0, 0, 0]], "CCOO")

aceticAcidSMILES = "CC(O)=O"

--cyclicMol :: ([[Int]], [Atom])
cyclicMol = ([[0, 1, 0, 0, 1, 0],
              [1, 0, 1, 0, 0, 0],
              [0, 1, 0, 1, 0, 0],
              [0, 0, 1, 0, 1, 0],
              [1, 0, 0, 1, 0, 1],
              [0, 0, 0, 0, 1, 0]] , "CCCNCO")

cmEL = adjToEdges (fst cyclicMol)
aaEL = adjToEdges (fst aceticAcid)
cmG  = toGraph cmEL
aaG  = toGraph aaEL

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

toGraph el = buildG (0, m) el
              where m = maximum $ nub (map fst el ++ map snd el)

-- build up a minmum spanning tree
prims :: [(Int, Int)] -> [(Int, Int)]
prims [] = []
prims el = sort $ prims' m [0] [] el
      where m = maximum $ concat [[fst x, snd x] | x <- el]

prims' m ns es el | length ns == m + 1 = es
                  | otherwise = prims' m ((fst e):ns) ((snd e):es) el
                      where e = findEdge ns el

findEdge      :: [Int] -> [(Int, Int)] -> (Int, (Int, Int))
findEdge ns [] = error "No edges found"
findEdge ns (e:es) | (elem (fst e) ns) && (notElem (snd e) ns) = ((snd e), e)
                   | (notElem (fst e) ns) && (elem (snd e) ns) = ((fst e), e)
                   | otherwise = findEdge ns es

cyclicEdges el = [e | e <- el, notElem e p]
      where p = prims el

bondCount el = [countBond e el | e <- p]
      where p = prims el

countBond e el = (c, e)
      where c = length [e' | e' <- el, e' == e]

contents i []     = []
contents i (a:as) | i == fst a    = snd a
                  | otherwise = contents i as

makeSmiles mol = smilesify mol 0

smilesify mol n = [key !! n] ++ ringChk n cyc bc ++ smileMore n (contents n al) mol
      where key = snd mol
            el = adjToEdges (fst mol)
            el' = simpleEdgeList el
            p = prims el
            al = edgeListToAdjList p
            bc  = bondCount el
            cyc = cyclicEdges el
            m = length key

smileMore p [] mol     = ""
smileMore p (v:vs) mol | null vs   = bondChk (p, v) bc ++ smilesify mol v
                       | otherwise =  "(" ++ bondChk (p, v) bc ++ smilesify mol v ++ ")" ++ smileMore p vs mol
                        where el = adjToEdges (fst mol)
                              bc  = bondCount el
                              cyc = cyclicEdges el

--bondChk :: (Int, Int) -> [(Int, (Int, Int))] -> [Char]
bondChk e []                  = ""
bondChk e (b:bs) | e == snd b = bondDict (fst b)
                 | otherwise  = bondChk e bs

--ringChk :: Int -> [(Int, Int)] -> [(Int, (Int, Int))] -> [Char]
ringChk n [] _               = ""
ringChk n cyc bc | elem n v  = concat
                              [show i ++ bondChk e bc |
                              (i, e) <- zip [0..] cyc, n == fst e || n == snd e]
                 | otherwise = ""
                 where v = nub (map fst cyc ++ map snd cyc)
