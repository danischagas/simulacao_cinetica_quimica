import numpy as np
import random as rd
from itertools import combinations

######################################################################################################################################################################
#                                                                                                                                                                    #
#                                                                     Classe                                                                                         #
#                                                                                                                                                                    #
######################################################################################################################################################################


class Particula:
    """Define physics of elastic colisao."""
    
    def __init__(self, massa, raio, posicao, velocidade):
        """Initialize a Particle object
        
        massa the massa of particle
        raio the raio of particle
        posicao the posicao vector of particle
        velocidade the velocidade vector of particle
        """
        self.massa = massa
        self.raio = raio
        self.tipo = 'atomo'
        
        # last posicao and velocidade
        self.posicao = np.array(posicao)
        self.velocidade = np.array(velocidade)
        
        # all posicao and velocities recorded during the simulation
        self.todas_posicoes = [np.copy(self.posicao)]
        self.todas_velocidades = [np.copy(self.velocidade)]
        self.todas_magnitudesv = [np.linalg.norm(np.copy(self.velocidade))]
        
    def prox_passo(self, passo):
        """Compute posicao of next passo."""
        self.posicao += passo * self.velocidade
        self.todas_posicoes.append(np.copy(self.posicao)) 
        self.todas_velocidades.append(np.copy(self.velocidade)) 
        self.todas_magnitudesv.append(np.linalg.norm(np.copy(self.velocidade))) 
    
    def checar_colisão(self, particle):
        """Check if there is a colisao with another particle."""
        
        r1, r2 = self.raio, particle.raio
        x1, x2 = self.posicao, particle.posicao
        di = x2-x1
        norm = np.linalg.norm(di)
        if norm-(r1+r2)*1.1 < 0:
            return True
        else:
            return False

    def realiza_colisao(self, particle, passo):
        """Compute velocidade after colisao with another particle."""
        m1, m2 = self.massa, particle.massa
        r1, r2 = self.raio, particle.raio
        v1, v2 = self.velocidade, particle.velocidade
        x1, x2 = self.posicao, particle.posicao
        di = x1-x2 # mudamos para ficar igual a fórmula
        norm = np.linalg.norm(di)
        if norm-(r1+r2)*1.1 < passo*abs(np.dot(v1-v2, di))/norm:
            self.velocidade = v1 - 2. * m2/(m1+m2) * np.dot(v1-v2, di) / (np.linalg.norm(di)**2.) * di
            particle.velocidade = v2 - 2. * m1/(m2+m1) * np.dot(v2-v1, (-di)) / (np.linalg.norm(di)**2.) * (-di)
            

    def colisao_paredes(self, passo, size):
        """Compute velocidade after hitting an edge.
        passo the computation passo
        size the medium size
        """
        r, v, x = self.raio, self.velocidade, self.posicao
        projx = passo*abs(np.dot(v,np.array([1.,0.])))
        projy = passo*abs(np.dot(v,np.array([0.,1.])))
        if abs(x[0])-r < projx or abs(size-x[0])-r < projx:
            self.velocidade[0] *= -1
        if abs(x[1])-r < projy or abs(size-x[1])-r < projy:
            self.velocidade[1] *= -1.
            
######################################################################################################################################################################
#                                                                                                                                                                    #
#                                                                    Funções                                                                                         #
#                                                                                                                                                                    #
######################################################################################################################################################################

def mudar_passo(lista_particulas, passo, size):
    """Solve a passo for every particle."""
    
    # Detect edge-hitting and colisao of every particle
    for i in range(len(lista_particulas)):
        lista_particulas[i].colisao_paredes(passo,size)
        for j in range(i+1,len(lista_particulas)):
                lista_particulas[i].realiza_colisao(lista_particulas[j],passo)    

                
    # Computa a posicao de todas particulas  
    for particula in lista_particulas:
        particula.prox_passo(passo)


def gerar_particulas(N, raio, massa, tamanho_caixa):
    """Generate N Particle objects in a random way in a list."""
    lista_particulas = []

    for i in range(N):
        
        magnitude_velocidade = np.random.rand(1)*25
        angulo_velocidade = np.random.rand(1)*2*np.pi
        v = np.append(magnitude_velocidade*np.cos(angulo_velocidade), magnitude_velocidade*np.sin(angulo_velocidade))
        
        colisao = True
        while(colisao == True):
            
            colisao = False
            pos = raio + np.random.rand(2)*(tamanho_caixa-2*raio) 
            nova_particula = Particula(massa, raio, pos, v)
            for j in range(len(lista_particulas)):

                colisao = nova_particula.checar_colisão( lista_particulas[j] )

                if colisao == True:
                    break

        lista_particulas.append(nova_particula)
    return lista_particulas


def simular_colisao(lista_particulas, probabilidade_reacao):
    for particula1, particula2 in combinations(lista_particulas, 2):
        if particula1.checar_colisão(particula2) and particula1.tipo == 'atomo' and particula2.tipo == 'atomo':
            valor_aleatorio = rd.random()

            if valor_aleatorio < probabilidade_reacao:
                nova_massa = particula1.massa + particula2.massa
                nova_posicao = (particula1.posicao + particula2.posicao) / 2 #centro de massa
                nova_velocidade = (particula1.massa * particula1.velocidade + particula2.massa * particula2.velocidade) / nova_massa
                nova_particula = Particula(
                    massa=nova_massa,
                    raio=(particula1.raio**3 + particula2.raio**3)**(1/3),  # Lei da conservação de volume
                    tipo='molecula',
                    posicao=nova_posicao,
                    velocidade=nova_velocidade
                ) 
                
                lista_particulas.remove(particula1)
                lista_particulas.remove(particula2)
                lista_particulas.append(nova_particula)
                return nova_particula
            else:
                #colisão elastica
                return realiza_colisao(particula1, particula2, passo)
    else:
        # Não houve colisão
        return None   
