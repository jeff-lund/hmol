{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}

import Data.Char
import Data.String
import Data.Void
import Text.Parsec
import Text.Parsec.Char
import Text.Parsec.Combinator
import Text.Parsec.String
import Text.Parsec.Token
import Data.Either
import Treedot

aceticAcid = "CC(O)=O"
caffeine   = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
cyclohexene = "C=1CCCCC=1"
long_chain = "CCCCC(CCO)C(C)OC"
halogen_start = "ClNO(C(C)O)CO(C)"
arm = "Cc1cccc1ON"
aliphtest = "CCO=CC=N=ClC=C"
newtest = "C(O(N(Cl(Br(IO)))))"

parseG rule text = extract $ parse rule "(source)" text
extract e = case e of Right e' -> e'
                      Left err -> error "Failed Parse"
-- organic subset
data AtomSymbol = B | C | N | O | S | P | F | Cl | Br | I
                  deriving (Show, Eq, Enum, Read)

data Atom = Aliphatic AtomSymbol | Aromatic AtomSymbol
            deriving (Eq, Read)

instance Show Atom where
    show (Aliphatic atom) = show atom
    show (Aromatic atom)  = map toLower (show atom)

data BondType = Single | Double | Triple | Quadruple | Arom
            deriving (Show, Eq, Enum, Read)

bondDict b | (b == '-') = Single
           | (b == '=') = Double
           | (b == '#') = Triple
           | (b == '$') = Quadruple

data Chain = Bond BondType
           | Atom Atom
           | RingBond (Maybe BondType) Int
            deriving (Eq, Read)
-- Comment this out and derive show above for old school dot graphs
instance Show Chain where
    show (Bond b) = show b
    show (Atom a) = show a
    show (RingBond b n) = case b of
                              Nothing -> "Ring " ++ show n
                              Just b' -> (show b') ++ " | Ring " ++ (show n)

data SMILESTree a = Node a [SMILESTree a]
            deriving (Show)

-- Treedot defs
instance LabeledTree (SMILESTree Chain)
    where label (Node x _) = (show x)
instance Tree (SMILESTree Chain)
    where subtrees (Node _ as) = as

bonds     = "-=#$"
atoms     = map show [B ..]
aromatics = concat [map toLower x | x <- (map show [B, C, N, O, S, P])]

twoAliphaticP :: Parser Chain
twoAliphaticP = do u <- upper
                   l <- oneOf "rl"
                   return (Atom (Aliphatic (read [u,l]::AtomSymbol)))

oneAliphaticP :: Parser Chain
oneAliphaticP = do u <- upper
                   return (Atom (Aliphatic (read [u]::AtomSymbol)))

aromaticP :: Parser Chain
aromaticP = do a <- oneOf aromatics
               return (Atom (Aromatic (read [toUpper a]::AtomSymbol)))

atomP :: Parser Chain
atomP = try twoAliphaticP <|> oneAliphaticP <|> aromaticP

bondP :: Parser Chain
bondP = do b <- oneOf bonds
           d <- optionMaybe (many1 digit)
           case d of
             Nothing -> return (Bond (bondDict b))
             Just d' -> return (RingBond (Just (bondDict b)) (read $ concat d::Int))

ringP :: Parser Chain
ringP = do d <- (many1 digit)
           return (RingBond Nothing (read d::Int))

parenP :: Parser (SMILESTree Chain)
parenP = do char '('
            r <- smilesP
            char ')'
            return r

smilesP :: Parser (SMILESTree Chain)
smilesP = do a <- try atomP <|> try ringP <|> bondP
             p <- optionMaybe (many parenP)
             r <- many smilesP
             case p of
               Nothing -> pure (Node a r)
               Just p' -> pure (Node a (p' ++ r))
