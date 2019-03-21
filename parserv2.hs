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

data BondType = Single
              | Double
              | Triple
              | Quadruple
              | RingBond (BondType, Ring)
            deriving (Show, Eq)

newtype Ring = Ring Int
    deriving (Show, Eq, Read)

{- instance Show Ring where
    show Ring (Nothing, i)  = show i
    show Ring ((Just b), i) = (show b) ++ " " ++ (show i) -}



newtype Chain = Atom' (Atom, (Maybe BondType), (Maybe Ring))

instance Show Chain where
  show (Atom' (a, Nothing, Nothing)) = show a
  show (Atom' (a, Just b, Nothing))  = (show b) ++ ", " ++ (show a)
  show (Atom' (a, Nothing, Just r))  = (show a) ++ ", " ++ (show r)
  show (Atom' (a, Just b, Just r))   = (show b) ++ ", " ++ (show a) ++ ", " ++ (show r)

data SMILESTree a = Node a [SMILESTree a]
  deriving (Show)

instance LabeledTree (SMILESTree Chain)
    where label (Node x _) = (show x)

instance Tree (SMILESTree Chain)
    where subtrees (Node _ as) = as

bonds     = "-=#$"
atoms     = map show [B ..]
aromatics = concat [map toLower x | x <- (map show [B, C, N, O, S, P])]
bondDict b | (b == '-') = Single
           | (b == '=') = Double
           | (b == '#') = Triple
           | (b == '$') = Quadruple

twoAliphaticP :: Parser Atom
twoAliphaticP =  do u <- upper
                    l <- oneOf "rl"
                    return (Aliphatic (read [u,l]::AtomSymbol))

oneAliphaticP :: Parser Atom
oneAliphaticP =  do u <- upper
                    return (Aliphatic (read [u]::AtomSymbol))

aromaticP :: Parser Atom
aromaticP =  do a <- oneOf aromatics
                return (Aromatic (read [toUpper a]::AtomSymbol))

atomP :: Parser Atom
atomP = try twoAliphaticP <|> try oneAliphaticP <|> aromaticP

bondP :: Parser BondType
bondP = do b <- oneOf bonds
           d <- optionMaybe (many1 digit)
           case d of
             Nothing -> return (bondDict b)
             Just d' -> return (RingBond (bondDict b, Ring (read $ concat d::Int)))

ringP :: Parser Ring
ringP = do d <- (many1 digit)
           return (Ring (read d::Int))

parenP :: Parser (SMILESTree Chain)
parenP = do char '('
            r <- smilesP
            char ')'
            return r

smilesP :: Parser (SMILESTree Chain)
smilesP = do b <- optionMaybe bondP
             a <- atomP
             r <- optionMaybe ringP
             p <- optionMaybe (many parenP)
             c <- many smilesP
             case p of
               Nothing -> pure (Node (Atom' (a, b, r)) c)
               Just p' -> pure (Node (Atom' (a, b, r)) (p' ++ c))


acidTree = parseG smilesP aceticAcid
