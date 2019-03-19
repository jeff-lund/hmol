import Data.Graph
import Data.List
import Data.Array

data Atom = H  | He | Li | Be | B  | C  | N  | O  | F  | Ne | Na | Mg | Al |
            Si | P  | S  | Cl | Ar | K  | Ca | Sc | Ti | V  | Cr | Mn | Fe |
            Co | Ni | Cu | Zn | Ga | Ge | As | Se | Br | Kr | Rb | Sr | Y  |
            Zr | Nb | Mo | Tc | Ru | Rh | Pd | Ag | Cd | In | Sn | Sb | Te |
            I  | Xe | Cs | Ba | Hf | Ta | W  | Re | Os | Ir | Pt | Au | Hg |
            Tl | Pb | Bi | Po | At | Rn | Fr | Ra | Rf | Db | Sg | Bh | Hs |
            Mt | Ds | Rg | Cn | Fl | Lv | La | Ce | Pr | Nd | Pm | Sm | Eu |
            Gd | Tb | Dy | Ho | Er | Tm | Yb | Lu | Ac | Th | Pa | U  | Np |
            Pu | Am | Cm | Bk | Cf | Es | Fm | Md | No | Lr | Star |
            B_| C_ | N_ | O_ | P_ | S_ | Se_ | As_ -- armoatics
            deriving (Show, Eq)

data Bond = Single | Double | Triple | Aromatic

organics = [B, C, N, O, S, P, F, Cl, Br, I, B_, C_, N_, O_, S_, P_, Star]

bondOrder     :: Atom -> Maybe [Int]
bondOrder a | a == C                = Just [4]
            | a == N                = Just [3, 5]
            | a == O                = Just [2]
            | a == P                = Just [3, 5]
            | a == S                = Just [2, 4, 6]
            | elem a [F, Br, Cl, I] = Just [1]
            | a == Star             = Just []
            | otherwise             = Nothing

--aceticAcid :: ([[Int]], [Atom])
aceticAcid = ([[0, 1, 2, 1],
               [1, 0, 0, 0],
               [2, 0, 0, 0],
               [1, 0, 0, 0]], [C, C, O, O])

aceticAcidSMILES = "CC(O)=O"

--cyclicMol :: ([[Int]], [Atom])
cyclicMol = ([[0, 1, 0, 0, 1, 0],
              [1, 0, 1, 0, 0, 0],
              [0, 1, 0, 1, 0, 0],
              [0, 0, 1, 0, 1, 0],
              [1, 0, 0, 1, 0, 1],
              [0, 0, 0, 0, 1, 0]] , [C, C, C, C, C, O])

cmEL = adjToEdges (fst cyclicMol)
aaEL = adjToEdges (fst aceticAcid)
cmG  = toGraph cmEL
aaG  = toGraph aaEL

adjToEdges []     = []
adjToEdges rs = convert 0 0 rs

adjToSimpleEdges [] = []
adjToSimpleEdges m  = convert' 0 m

convert n k [] = []
convert n k (r:rs) = concat [multiEdge c (n, ind) | (c, ind) <- zip r [0..], c > k, ind > n]
                     ++ convert (n+1) k rs


convert' n [] = []
convert' n (r:rs) = [(n, ind) | (c, ind) <- zip r [0..], c > 0, ind > n]
                     ++ convert' (n+1) rs

multiEdge 0 e = []
multiEdge n e =  e:(multiEdge (n-1) e)

toGraph edgeList = buildG (0, m) edgeList
              where m = maximum $ concat [[fst x, snd x] | x <- edgeList]

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
